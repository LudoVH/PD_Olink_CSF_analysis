
#load required libraries
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(stats)
library(reshape2)
library(tidyverse)
library(ggbiplot) 
library(readxl)
library(umap)
library(gtsummary)
library(labelled)
library(plotly)
library(table1)
library(writexl)
options(ggrepel.max.overlaps = Inf)

############ OLINK DATA ###############

# read NPX data into R 
PD_Olink_CSF_raw <- read_NPX("~/Documents/R/Proteomics/Olink/CSF/Olink_proteomics_CSF/data/P230168_bridged_NPX_data_adj_factor_intensity.csv") 
PD_Olink_CSF_raw <- PD_Olink_CSF_raw %>%
        mutate(SampleID = str_replace_all(SampleID, "PDS_PA_010_2", "PDS_PA_010_1")) # fixes a mistype 

# summary table of raw data
var_label(PD_Olink_CSF_raw$Assay_Warning) <- "Assay Warning"
var_label(PD_Olink_CSF_raw$QC_Warning) <- "QC Warning"
PD_Olink_CSF_raw_summary_table <- PD_Olink_CSF_raw %>%
        select(Sample_Type, Assay_Warning, QC_Warning)  %>%
        tbl_summary (missing_text = "Missing",
                     statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"))
print(PD_Olink_CSF_raw_summary_table)

# remove sample and assay controls 
PD_Olink_CSF_raw_without_controls <- PD_Olink_CSF_raw %>% 
        mutate(SampleID = str_replace_all(SampleID, "_", "/")) %>%
        filter(!str_detect(Sample_Type, "CONTROL")) %>%    # remove sample controls 
        filter(!str_detect(Assay, "control")) # remove assay controls
        
# create list of longitudinal participant IDs tested
Longitudinal_PD_samples <- c("PDS/JR/078/2", "PDS/JR/327/2", "PDS/JR/473/2")

# UMAP of raw data (Without sample & assay controls)
olink_umap_plot(PD_Olink_CSF_raw_without_controls, color_g = "QC_Warning", byPanel = FALSE)

############ MANIFEST DATA ###############

# import PD CSF Olink plating manifest variables
PD_Manifest_variables <- read_excel("~/Documents/R/Proteomics/Olink/CSF/CSF_Olink_manifest_FINAL.xlsx") %>%
        rename(dob = "Date of Birth", csf_date = "Biobanked CSF date", age_at_sample = "Age at biosampling", thaws = "Thaws (inclusive)") %>%
        select(SampleID, SubjectID, gender, cohort, diagnosis, csf_date, age_at_sample, thaws, plate) %>%
        mutate(csf_date = as.numeric(csf_date)) %>%
        mutate(csf_date = as.Date(csf_date, origin = "1899-12-30")) %>%
        mutate(date_of_testing = as.Date("2023-07-01")) %>%
        mutate(age_of_sample = round(time_length((date_of_testing - csf_date), unit = "years"))) %>%
        mutate(gender = as.character(gender)) %>%
        mutate(SampleID = str_replace_all(SampleID, "_", "/")) %>%
        mutate(disease_group = ifelse(cohort == "PD" & diagnosis == "HC", "HC-PD", diagnosis)) %>% 
        mutate(disease_group = ifelse(cohort == "ALS" & disease_group == "HC", "HC-ALS", disease_group)) %>%
        mutate(disease_state = ifelse(diagnosis == "HC", "control", "disease")) %>%
        mutate(thaws = as.character(thaws)) %>%
        mutate(LongitudinalSampleNo = ifelse(SampleID %in% Longitudinal_PD_samples, "2", "1")) 

# Relevant OPDC clinical data imported of the tested PD samples 
Relevant_OPDC_clinical_data <- read.csv("~/Documents/R/OPDC/OPDC Clinical data/OPDCDiscoveryCohort_DATA_2024-01-29_1345.csv") %>%
        select(subjid, dob, gender, visit_date, x_1_onest_date, x_2_diag_date, x_1_1_rbd_symptom_date, x_1_2_rbd_diag_date, ak_1_has_pd, ak_2_meets_criteria, converted_to_pd, converted_to_pd_r, converted_to_other, converted_to_other_diag,converted_subjid_2, alternate_diag, alternate_diag_1, alternate_diag_specify) %>%
        rename(SubjectID = subjid, date_onset_PD = x_1_onest_date, date_diagnosis_PD = x_2_diag_date, likelihood_PD = ak_1_has_pd, meets_UKBB = ak_2_meets_criteria, old_ID = converted_subjid_2, date_onset_RBD = x_1_1_rbd_symptom_date, date_diagnosis_RBD = x_1_2_rbd_diag_date) %>%
        mutate(dob = as.Date(dob), visit_date = as.Date(visit_date), date_onset_PD = as.Date(date_onset_PD), date_diagnosis_PD = as.Date(date_diagnosis_PD), date_onset_RBD = as.Date(date_onset_RBD), date_diagnosis_RBD = as.Date(date_diagnosis_RBD)) %>%
        mutate(gender = replace(gender, gender == 0, "male"), gender = replace(gender, gender == 1, "female")) %>%
        filter(SubjectID %in% PD_Manifest_variables$SubjectID)

# Date of disease diagnosis
Date_of_diagnosis <- bind_rows(Relevant_OPDC_clinical_data[!is.na(Relevant_OPDC_clinical_data$date_diagnosis_PD), ], Relevant_OPDC_clinical_data[!is.na(Relevant_OPDC_clinical_data$date_diagnosis_RBD), ]) %>%
                mutate(date_diagnosis = as.Date(ifelse(is.na(date_diagnosis_PD), date_diagnosis_RBD, date_diagnosis_PD))) %>%
                mutate(date_onset = as.Date(ifelse(is.na(date_onset_PD), date_onset_RBD, date_onset_PD))) %>%
                select(SubjectID, date_diagnosis, date_onset)

# add date of disease diagnosis to manifest
PD_Manifest_variables <- left_join(PD_Manifest_variables, Date_of_diagnosis, by = "SubjectID")

# import PD cluster assignment and add to manifest variables 
OPDC_PD_clusters <- read_excel("~/Documents/R/OPDC/Discovery_PD_clusters.xlsx")                              
OPDC_PD_clusters <- OPDC_PD_clusters %>% rename (SubjectID = SubjID, cluster = clusters)       
PD_Manifest_variables <- left_join(PD_Manifest_variables, OPDC_PD_clusters, by = "SubjectID")

# summary table of manifest
PD_Manifest_variables_table <- PD_Manifest_variables %>%
        tbl_summary(by = cohort, 
                    include = c(gender, diagnosis, LongitudinalSampleNo), 
                    label = c(diagnosis ~ "Diagnosis", gender ~ "Gender", LongitudinalSampleNo ~ "Longitudinal Sample No."),
                    statistic = c(diagnosis ~ "{n} ({p}%)", LongitudinalSampleNo ~ "{n} ({p}%)")) %>%
        add_overall() %>%
        modify_header(label ~ "**Cohort**") %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        as_gt() %>%
        gt::tab_source_note(gt::md("PD group includes 3 longitudinal samples.")) %>%
        gt::tab_source_note(gt::md("*ALS = amyotrophic lateral scleoris, DC = disease control, FDR = first-degree relative, HC = healthy control*")) %>%
        gt::tab_source_note(gt::md("*PD = Parkinson's disease, PLS = primary lateral sclerosis, RBD = REM-sleep behaviour disorder*"))
print(PD_Manifest_variables_table)

test <- PD_Manifest_variables %>%
        filter(diagnosis == "PD" | diagnosis == "RBD")

setdiff(Date_of_diagnosis$SubjectID, test$SubjectID)

################### RAW OLINK + MANIFEST DATA ############################

# Raw Olink results with variables (and HC warnings labeled)
PD_Olink_CSF_data_raw_variables <- PD_Olink_CSF_raw_without_controls %>% 
        select(-Sample_Type) %>%
        relocate(SampleID, .before = Index) %>%
        left_join(PD_Manifest_variables, by = "SampleID") %>%
        group_by(SampleID) %>%
        mutate(QC_WARN_percentage = sum(QC_Warning == "WARN")/sum(QC_Warning != "")*100) %>%
        mutate(QC_rate = ifelse(QC_WARN_percentage >=1, "High QC WARN %", "Low QC WARN %")) %>%
        ungroup()

# summary table of the raw data with variables 
var_label(PD_Olink_CSF_data_raw_variables$Assay_Warning) <- "Assay Warning"
var_label(PD_Olink_CSF_data_raw_variables$QC_Warning) <- "QC Warning"
PD_Olink_CSF_data_raw_variables_summary_table <- PD_Olink_CSF_data_raw_variables %>%
        select(diagnosis, cohort, Assay_Warning, QC_Warning)  %>%
        tbl_summary (by = diagnosis, missing_text = "Missing",
                     statistic = list(cohort = "{n}", Assay_Warning = "{n} ({p}%)", QC_Warning = "{n} ({p}%)"))
print(PD_Olink_CSF_data_raw_variables_summary_table)

# UMAP of raw data with variable
olink_umap_plot(PD_Olink_CSF_data_raw_variables, color_g = "diagnosis", byPanel = FALSE)
olink_umap_plot(PD_Olink_CSF_data_raw_variables, color_g = "cohort")

################### EXPLORATION OF SAMPLES WITH HIGH QC WARNING PERCENTAGE ############################

# UMAP of raw data - QC warning 
olink_umap_plot(PD_Olink_CSF_data_raw_variables, color_g = "QC_Warning", byPanel = FALSE)
olink_qc_plot(PD_Olink_CSF_data_raw_variables, color_g = "QC_Warning")

QC_Warnings_PD <- PD_Olink_CSF_data_raw_variables %>%
        group_by(SampleID) %>%
        distinct(SampleID,.keep_all = TRUE) %>%
        select(SampleID, SubjectID, cohort, diagnosis, gender, QC_WARN_percentage, age_of_sample, LongitudinalSampleNo, cluster, thaws, QC_rate) 

# density plot of QC_WARN_percentage 
plot(density(QC_Warnings_PD$QC_WARN_percentage, bw = 0.1), xlab = "QC WARN percentage", main = "Density plot of QC WARN percentage")

# Selects out the SampleIDs with a QC_WARN_percentage of 1% or greater
High_QC_Warnings_PD <- QC_Warnings_PD %>% 
        filter(QC_WARN_percentage >= 1)

# Summary table of high QC warnings 
High_QC_Warnings_PD_table <- High_QC_Warnings_PD %>%
        tbl_summary(by = cohort,
                    include = c(diagnosis),
                    label = list(diagnosis ~ "Diagnosis"),
                    statistic = list(cohort = "{n}")) %>%
        add_overall()
print(High_QC_Warnings_PD_table)

write_xlsx(High_QC_Warnings_PD,"~/Documents/R/Proteomics/Olink/CSF/High_QC_outliers.xlsx") # save list to excel file


######################## DATA CLEANED ######################

# Remove all readouts with an Olink-generated QC Warning that is "EXCLUDED" or "WARN"
PD_Olink_CSF_data_clean <- PD_Olink_CSF_data_raw_variables %>% 
        filter(!str_detect(QC_Warning, "EXCLUDED|WARN")) %>% # remove results with sample QC warnings 
        filter(!str_detect(Assay_Warning, "EXCLUDED|WARN")) # remove assay QC warnings 

olink_umap_plot(PD_Olink_CSF_data_clean)

# UMAP of cleaned data
olink_umap_plot(PD_Olink_CSF_data_clean, color_g = "QC_rate", byPanel = FALSE)
olink_umap_plot(PD_Olink_CSF_data_clean, color_g = "diagnosis", byPanel = FALSE)
olink_dist_plot(PD_Olink_CSF_data_clean, color_g = "QC_rate")

PD_Olink_CSF_data_clean %>%
        mutate(Panel = "Olink") %>%
        olink_dist_plot(color_g = "QC_rate")

PD_Olink_CSF_data_clean_without_bad_samples <- PD_Olink_CSF_data_clean %>%
        filter(!QC_rate == "High QC WARN %") 
olink_umap_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "cohort")


############################## EXPLORING NIKOLETA'S OUTLIERS ################################

# import outliers from Nikoleta 
Nikoleta_outliers <- c("PDS-AH-010_1", "PDS-AH-074_1", "PDS-AH-083_1", "PDS-AH-102_1", "PDS-JR-007_1", "PDS-JR-015_1", "PDS-JR-022_1", "PDS-JR-023_1", "PDS-JR-024_1",
                       "PDS-JR-029_1", "PDS-JR-034_1", "PDS-JR-092_1", "PDS-JR-147_1", "PDS-JR-189_1", "PDS-JR-230_1", "PDS-JR-261_1", "PDS-JR-334_1", "PDS-JR-340_1",
                       "PDS-JR-341_1", "PDS-JR-354_1", "PDS-JR-423_1", "PDS-JR-424_1", "PDS-MK-014_1", "PDS-MK-017_1", "PDS-MK-020_1", "PDS-MK-056_1", "PDS-RH-004_1",
                       "PDS-RH-007_1", "PDS-RH-024_1", "PDS-WP-036_1", "PDS-WP-047_1") %>%
        str_replace_all("-", "/") %>%
        str_replace_all("_", "/")

# data of Nikoleta's outliers 
PD_Nikoleta_outliers_data <- PD_Olink_CSF_data_clean %>%
        filter(SampleID %in% Nikoleta_outliers) %>%
        group_by(SampleID) %>%
        distinct(SampleID, .keep_all = TRUE)

# summary table 1 of Nikoleta's outliers 
PD_Nikoleta_outliers_table1 <- PD_Nikoleta_outliers_data %>%
        tbl_summary(include = c(cohort, diagnosis, gender, age_of_sample),
                    label = list(cohort ~ "Cohort", diagnosis ~ "Diagnosis", gender ~ "Gender", age_of_sample ~ "Age of sample")) %>%
        modify_header(label ~ "**Variable**") 
print(PD_Nikoleta_outliers_table1)

# summary table 2 of Nikoleta's outliers 
PD_outliers_table2 <- PD_Nikoleta_outliers_data %>%
        tbl_summary(include = c(LongitudinalSampleNo, cluster, thaws, QC_rate),
                    label = list(LongitudinalSampleNo ~ "Longitudinal Sample No.", thaws ~ "No. of freeze/thaw cycles", cluster ~ "PD cluster", QC_rate ~ "Rate of QC Warnings")) %>%
        modify_header(label ~ "**Variable**") 
print(PD_outliers_table2)


# Plots with Nikoleta's outliers labelled
PD_Olink_CSF_data_clean_with_outliers <- PD_Olink_CSF_data_clean %>%
        mutate(outlier = ifelse(SampleID %in% Nikoleta_outliers, "yes(1)", "no(1)"))

olink_umap_plot(PD_Olink_CSF_data_clean_with_outliers, color_g = "outlier", byPanel = FALSE)
olink_qc_plot(PD_Olink_CSF_data_clean_with_outliers, color_g = "outlier")
olink_dist_plot(PD_Olink_CSF_data_clean_with_outliers, color_g = "outlier")
olink_heatmap_plot(PD_Olink_CSF_data_clean_with_outliers)

PD_Olink_CSF_data_clean_with_outliers %>% 
        mutate(Panel = "Olink") %>%
        olink_dist_plot(color_g = "outlier")

# remove Nikoleta's outliers
PD_Olink_CSF_data_clean_without_Nikoleta_outliers <- PD_Olink_CSF_data_clean %>%
        filter(!SampleID %in% Nikoleta_outliers)

write_xlsx(PD_Nikoleta_outliers_data,"~/Documents/R/Proteomics/Olink/CSF/Nikoleta_outliers.xlsx") # save list to excel file


######################################## FINAL DATA SET #################################################

# Cleaned data with bad samples removed 
PD_Olink_CSF_data_clean_without_bad_samples <- PD_Olink_CSF_data_clean %>%
        filter(!SampleID %in% Nikoleta_outliers)

# write to csv file 
export(PD_Olink_CSF_data_clean_without_bad_samples, "PD_Olink_CSF_data_clean_without_bad_samples.csv")


# Scatterplot
scatterplot_QC <- olink_qc_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "diagnosis")
scatterplot_outliers <- scatterplot_QC$data %>% filter(Outlier == 1) %>% select(SampleID, Panel, diagnosis, IQR, sample_median, Outlier)

# summary table of QC outliers
scatterplot_outliers_table <- scatterplot_outliers %>%
        tbl_summary (by = SampleID, 
                     include = c(Panel, diagnosis),
                    statistic = Panel ~ "{n}")
print(scatterplot_outliers_table)


# Summary table of final data set 
PD_Olink_CSF_data_final_table <- PD_Olink_CSF_data_clean_without_bad_samples %>%
        group_by(SampleID) %>%
        distinct(SampleID, .keep_all = TRUE) %>%
        tbl_summary(by = cohort, 
                    include = c(diagnosis, gender)) %>%
        add_overall() %>%
        modify_header(label ~ "**Cohort**") %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        as_gt() %>%
        gt::tab_source_note(gt::md("PD group includes 3 longitudinal samples.")) %>%
        gt::tab_source_note(gt::md("*ALS = amyotrophic lateral scleoris, DC = disease control, FDR = first-degree relative, HC = healthy control*")) %>%
        gt::tab_source_note(gt::md("*PD = Parkinson's disease, PLS = primary lateral sclerosis, RBD = REM-sleep behaviour disorder*"))
print(PD_Olink_CSF_data_final_table)


# Final UMAP and PCA of data
olink_umap_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "diagnosis", byPanel = FALSE)
olink_umap_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "diagnosis", byPanel = TRUE)
olink_pca_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "diagnosis", byPanel = FALSE)
olink_pca_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "diagnosis", byPanel = TRUE)
olink_umap_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "cohort", byPanel = FALSE)
olink_pca_plot(PD_Olink_CSF_data_clean_without_bad_samples, color_g = "cohort", byPanel = FALSE)


######################## COMPARING THE HC DATA ######################

# HC-PD vs HC-ALS 
PD_Olink_CSF_HC_only <- PD_Olink_CSF_data_clean_without_bad_samples %>% 
        filter(diagnosis == "HC")
olink_umap_plot(PD_Olink_CSF_HC_only, color_g = "disease_group", byPanel = FALSE) 
olink_pca_plot(PD_Olink_CSF_HC_only, color_g = "disease_group", byPanel = FALSE) 

olink_qc_plot(PD_Olink_CSF_HC_only, color_g = "disease_group")




# perform t-test
ttest_results_HCvsHC <- olink_ttest(df = PD_Olink_CSF_HC_only,
                                    variable = 'disease_group')
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_HCvsHC,
                   x_lab = 'HC-PD vs HC-ALS',
                   olinkid_list = FALSE)

                
######################## PD COHORT DATA ######################

# Final data set of PD data only         
PD_Olink_CSF_data_PD_cohort <- PD_Olink_CSF_data_clean_without_bad_samples %>%
        filter(cohort == "PD") 

# Summary table of PD data set 
PD_Olink_CSF_data_PD_cohort_table <- PD_Olink_CSF_data_PD_cohort %>%
        group_by(SampleID) %>%
        distinct(SampleID, .keep_all = TRUE) %>%
        tbl_summary(by = diagnosis, 
                    include = c(gender)) %>%
        add_overall() %>%
        modify_header(label ~ "**Cohort**") %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        as_gt() %>%
        gt::tab_source_note(gt::md("PD group includes 3 longitudinal samples.")) %>%
        gt::tab_source_note(gt::md("*ALS = amyotrophic lateral scleoris, DC = disease control, FDR = first-degree relative, HC = healthy control*")) %>%
        gt::tab_source_note(gt::md("*PD = Parkinson's disease, PLS = primary lateral sclerosis, RBD = REM-sleep behaviour disorder*"))
print(PD_Olink_CSF_data_PD_cohort_table)


# Final data set of PD data, excluding 3 longitudinal samples 
PD_Olink_CSF_data_PD_cohort_single_visit <- PD_Olink_CSF_data_clean_without_bad_samples %>%
        filter(cohort == "PD") %>%
        filter(LongitudinalSampleNo == 1)
        
# UMAP of PD cohort (excl longitudinal visits)
olink_umap_plot(PD_Olink_CSF_data_PD_cohort_single_visit, color_g = "diagnosis", byPanel = FALSE)
olink_pca_plot(PD_Olink_CSF_data_PD_cohort_single_visit, color_g = "diagnosis", byPanel = FALSE)





