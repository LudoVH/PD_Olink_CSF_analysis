
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

############ OLINK DATA ###############

# read NPX data into R 
PD_Olink_CSF_raw <- read_NPX("~/Documents/R/Proteomics/Olink/CSF/Olink_proteomics_CSF/data/P230168_bridged_NPX_data_adj_factor_intensity.csv") 
PD_Olink_CSF_raw <- PD_Olink_CSF_raw %>% 
        mutate(SampleID = str_replace_all(SampleID, "_", "/")) %>%
        mutate(SampleID = str_replace_all(SampleID, "PDS/PA/010/2", "PDS/PA/010/1"))
        

# summary table of data, excluding plate controls 
var_label(PD_Olink_CSF_raw$Assay_Warning) <- "Assay Warning"
var_label(PD_Olink_CSF_raw$QC_Warning) <- "QC Warning"
PD_Olink_CSF_raw_summary_table <- PD_Olink_CSF_raw %>%
        filter(!str_detect(Sample_Type, "CONTROL")) %>%  #exclude plate controls
        select(SampleID, Assay_Warning, QC_Warning)  %>%
        tbl_summary (missing_text = "Missing",
                     statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"))
print(PD_Olink_CSF_raw_summary_table)


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
        mutate(diagnosis = ifelse(diagnosis == "DC", "HC", diagnosis)) %>%
        mutate(disease_group = ifelse(cohort == "PD" & diagnosis == "HC", "HC-PD", diagnosis)) %>% 
        mutate(disease_group = ifelse(cohort == "ALS" & disease_group == "HC", "HC-ALS", disease_group)) %>%
        mutate(disease_state = ifelse(diagnosis == "HC", "control", "disease"))

# create list of participant IDs tested
PD_Manifest_ID_list <- unique(PD_Manifest_variables$SubjectID)
Duplicate_PD_samples <- c("PDS/JR/078", "PDS/JR/327", "PDS/JR/473")

# summary table of manifest
PD_Manifest_variables_table <- PD_Manifest_variables %>%
        tbl_summary(by = cohort, include = diagnosis, statistic = diagnosis ~ "{n}")
print(PD_Manifest_variables_table)

# import PD cluster assignment and add to manifest variables 
OPDC_PD_clusters <- read_excel("~/Documents/R/OPDC/Discovery_PD_clusters.xlsx")                              
OPDC_PD_clusters <- OPDC_PD_clusters %>% rename (SubjectID = SubjID, cluster = clusters)       
PD_Manifest_variables <- left_join(PD_Manifest_variables, OPDC_PD_clusters, by = "SubjectID")


################### RAW OLINK + MANIFEST DATA ############################

# Raw Olink results with variables (and HC warnings labeled)
PD_Olink_CSF_data_raw_variables <- PD_Olink_CSF_raw %>% 
        filter(!str_detect(Sample_Type, "CONTROL")) %>%  #exclude plate controls
        select(-Sample_Type) %>%
        relocate(SampleID, .before = Index) %>%
        left_join(PD_Manifest_variables, by = "SampleID") %>%
        group_by(SampleID) %>%
        mutate(QC_WARN_percentage = sum(QC_Warning == "WARN")/sum(QC_Warning != "")*100) %>%
        mutate(outlier = ifelse(QC_WARN_percentage >=1, "High QC WARN%", "Low QC WARN %")) %>%
        ungroup()

# summary table of the raw data with variables 
var_label(PD_Olink_CSF_data_raw_variables$Assay_Warning) <- "Assay Warning"
var_label(PD_Olink_CSF_data_raw_variables$QC_Warning) <- "QC Warning"
PD_Olink_CSF_data_raw_variables_summary_table <- PD_Olink_CSF_data_raw_variables %>%
        select(diagnosis, cohort, Assay_Warning, QC_Warning)  %>%
        tbl_summary (by = diagnosis, missing_text = "Missing",
                     statistic = list(cohort = "{n}", Assay_Warning = "{n} ({p}%)", QC_Warning = "{n} ({p}%)"))
print(PD_Olink_CSF_data_raw_variables_summary_table)


######################## DATA CLEANED ##################

# select all valid results 
PD_Olink_CSF_data_clean <- PD_Olink_CSF_data_raw_variables %>% 
        filter(!str_detect(QC_Warning, "EXCLUDED|WARN")) %>%
        filter(!str_detect(Assay_Warning, "EXCLUDED|WARN"))

# cleaned data with problem samples removed 
PD_Olink_CSF_without_outliers <- PD_Olink_CSF_data_clean %>%
        filter(QC_WARN_percentage < 1)

################### EXPLORATION OF SAMPLES WITH HIGH QC WARNING PERCENTAGE ############################

# UMAP of raw data - Diagnosis
test2 <- PD_Olink_CSF_data_raw_variables %>% 
        filter(str_detect(diagnosis, "PD|RBD|HC"))
olink_umap_plot(test2, color_g = "diagnosis", byPanel = FALSE)


# Selects out the SampleIDs with a QC_WARN_percentage of 1% or greater
High_QC_Warnings_PD <- PD_Olink_CSF_data_raw_variables %>% 
        filter(QC_WARN_percentage >= 1) %>%
        group_by(SampleID) %>%
        distinct(SampleID,.keep_all = TRUE) %>%
        select(SampleID, SubjectID, cohort, diagnosis, QC_WARN_percentage)

High_QC_Warnings_PD_table <- High_QC_Warnings_PD %>%
        tbl_summary(by = diagnosis,
                    include = c(cohort, diagnosis),
                    statistic = list(cohort = "{n}"))
print(High_QC_Warnings_PD_table)
        
# Creates character vector with all the High QC warning IDS
Problem_samples <- unique(High_QC_Warnings_PD$SampleID)




# UMAP of raw data - High QC rate
UMAP_Dx <- olink_umap_plot(High_QC_Warnings_PD, color_g = , byPanel = FALSE)


# UMAP of raw data by plate
UMAP_raw_panel <- olink_umap_plot(PD_Olink_CSF_data_raw_variables, color_g = "plate", byPanel = FALSE)

################### EXPLORATION OF OUTLIERS ############################

# import PD outliers from Nikoleta 
PD_outliers <- c("PDS-AH-010_1", "PDS-AH-074_1", "PDS-AH-083_1", "PDS-AH-102_1", "PDS-JR-007_1", "PDS-JR-015_1", "PDS-JR-022_1", "PDS-JR-023_1", "PDS-JR-024_1",
                 "PDS-JR-029_1", "PDS-JR-034_1", "PDS-JR-092_1", "PDS-JR-147_1", "PDS-JR-189_1", "PDS-JR-230_1", "PDS-JR-261_1", "PDS-JR-334_1", "PDS-JR-340_1",
                 "PDS-JR-341_1", "PDS-JR-354_1", "PDS-JR-423_1", "PDS-JR-424_1", "PDS-MK-014_1", "PDS-MK-017_1", "PDS-MK-020_1", "PDS-MK-056_1", "PDS-RH-004_1",
                 "PDS-RH-007_1", "PDS-RH-024_1", "PDS-WP-036_1", "PDS-WP-047_1") %>%
        str_replace_all("-", "/") %>%
        str_replace_all("_", "/")

# create table with Nikoleta's PD outliers 
PD_outliers_table <- PD_Olink_CSF_raw %>%
        mutate(SampleID = str_replace_all(SampleID, "_", "/")) %>%
        filter(SampleID %in% PD_outliers) %>%
        group_by(SampleID) %>%
        mutate(QC_WARN_percentage = sum(QC_Warning == "WARN")/sum(QC_Warning != "")*100) %>%
        mutate(QC_EXCLUDED_percentage = sum(QC_Warning == "EXCLUDED")/sum(QC_Warning != "")*100) %>%
        select(SampleID, QC_WARN_percentage) %>%
        distinct(SampleID, .keep_all = TRUE) 





########## EXPLORATION OF CLEANED DATA ######

# Diagnosis WITH outliers
UMAP_Dx <- olink_umap_plot(PD_Olink_CSF_data_clean, color_g = "diagnosis", byPanel = FALSE)
UMAP_Dx_panel <- olink_umap_plot(PD_Olink_CSF_data_clean, color_g = "diagnosis", byPanel = TRUE)

# PD + RBD vs HC 
test <- PD_Olink_CSF_data_clean %>%
        filter(str_detect(diagnosis, "PD|RBD|HC"))
        
# perform t-test
ttest_results_clean <- olink_ttest(df = test,
                                 variable = 'PD_prodromal_PD')
# select names of proteins to show
top_10_name_clean <- ttest_results_clean %>%
        slice_head(n = 10) %>%
        pull(OlinkID)
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_clean,
                   x_lab = 'Cohort',
                   olinkid_list = top_10_name_clean)


# High QC warning WITH outliers
UMAP_QC <- olink_umap_plot(PD_Olink_CSF_data_clean, color_g = "QC_WARN_rating", byPanel = FALSE)
UMAP_QC_panel <- olink_umap_plot(PD_Olink_CSF_data_clean, color_g = "QC_WARN_rating", byPanel = TRUE)


# UMAP of cleaned data WITHOUT outliers by plate
UMAP_raw_panel <- olink_umap_plot(PD_Olink_CSF_without_outliers, color_g = "plate", byPanel = FALSE)


# Diagnosis WITHOUT outliers
UMAP_DxWO <- olink_umap_plot(PD_Olink_CSF_without_outliers, color_g = "diagnosis", byPanel = FALSE)
UMAP_DxWO_panel <- olink_umap_plot(PD_Olink_CSF_without_outliers, color_g = "diagnosis", byPanel = TRUE)


# Cohort
UMAP_cohort_WO <- olink_umap_plot(PD_Olink_CSF_without_outliers, color_g = "cohort", byPanel = FALSE)
UMAP_cohort_WO_panel <- olink_umap_plot(PD_Olink_CSF_without_outliers, color_g = "cohort", byPanel = TRUE)

# PD Clusters
PD_Olink_CSF_clusters <- PD_Olink_CSF_without_outliers %>%
        filter(!str_detect(diagnosis, "RBD|HC")) %>%
        mutate(cluster = ifelse(is.na(cluster), "No cluster", cluster)) %>%
        filter(!cluster == "No cluster")
UMAP_clusters <- olink_umap_plot(PD_Olink_CSF_clusters, color_g = "cluster", byPanel = FALSE)
length(unique(PD_Olink_CSF_clusters$SampleID))

PCA_clusters <- olink_pca_plot(PD_Olink_CSF_clusters, color_g = "cluster", byPanel = TRUE)


# PD vs RBD (no HC)
PD_Olink_CSF_noHC <- PD_Olink_CSF_without_outliers %>%
        filter(str_detect(diagnosis, "PD|RBD")) %>%
olink_umap_plot(PD_Olink_CSF_noHC, color_g = "diagnosis", byPanel = FALSE)
olink_umap_plot(PD_Olink_CSF_noHC, color_g = "diagnosis", byPanel = TRUE)

olink_pca_plot(PD_Olink_CSF_noHC, color_g = "diagnosis", byPanel = FALSE)


# PD + RBD vs HC
PD_Olink_CSF_PD_RBDvHC <- PD_Olink_CSF_without_outliers %>%
        filter(str_detect(diagnosis,"PD|RBD|HC"))
olink_umap_plot(PD_Olink_CSF_PD_RBDvHC, color_g = "PD_prodromal_PD", byPanel = FALSE)
olink_umap_plot(PD_Olink_CSF_PD_RBDvHC, color_g = "PD_prodromal_PD", byPanel = TRUE)
olink_pca_plot(PD_Olink_CSF_PD_RBDvHC, color_g = "PD_prodromal_PD", byPanel = FALSE)

# perform t-test
ttest_results_PDRBDvHC <- olink_ttest(df = PD_Olink_CSF_PD_RBDvHC,
                             variable = 'PD_prodromal_PD')
# select names of proteins to show
top_10_name <- ttest_results_PDRBDvHC %>%
        slice_head(n = 10) %>%
        pull(OlinkID)
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_PDRBDvHC,
                   x_lab = 'Disease status',
                   olinkid_list = top_10_name)

#### Aromatic-L-amino-acid decarboxylase		DDC	Cardiometabolic
#### Glutamate decarboxylase 2		GAD2	Inflammation_II
#### Low affinity immunoglobulin epsilon Fc receptor		FCER2	Neurology


# PD v HC
PD_Olink_CSF_PDvHC <- PD_Olink_CSF_without_outliers %>%
        filter(str_detect(diagnosis,"PD|HC"))

# perform t-test
ttest_results_PDvHC <- olink_ttest(df = PD_Olink_CSF_PDvHC,
                             variable = 'disease_state')
# select names of proteins to show
top_10_name <- ttest_results_PDvHC %>%
        slice_head(n = 10) %>%
        pull(OlinkID)
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_PDvHC,
                   x_lab = 'Disease status',
                   olinkid_list = top_10_name)


# RBD v HC
PD_Olink_CSF_RBDvHC <- PD_Olink_CSF_without_outliers %>%
        filter(str_detect(diagnosis,"RBD|HC"))

# perform t-test
ttest_results_RBDvHC <- olink_ttest(df = PD_Olink_CSF_RBDvHC,
                                   variable = 'disease_state')
# select names of proteins to show
top_10_name <- ttest_results_RBDvHC %>%
        slice_head(n = 10) %>%
        pull(OlinkID)
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_RBDvHC,
                   x_lab = 'Disease status',
                   olinkid_list = top_10_name)


########### HEALTHY CONTROLS #####################


# HC-PD vs HC-ALS 
PD_Olink_CSF_HC_only <- PD_Olink_CSF_data_clean %>%
        filter(str_detect(diagnosis, "HC")) 
olink_umap_plot(PD_Olink_CSF_HC_only, color_g = "disease_group", byPanel = FALSE)
olink_pca_plot(PD_Olink_CSF_HC_only, color_g = "disease_group", byPanel = FALSE)

length(unique(PD_Olink_CSF_HC_only$SampleID)) # number of HC samples in the data set is 38 (21 PD HC, 17 ALS HC)
PD_Olink_CSF_HC_only %>%
        group_by(diagnosis) %>%
        mutate(unique_types = n_distinct(SampleID)) %>%
        select(unique_types)

PD_HC_list_PD <- PD_Olink_CSF_HC_only %>% filter(str_detect(diagnosis, "HC-PD"))
PD_HC_list_PD <- unique(PD_HC_list_PD$SampleID)
PD_HC_list_ALS <- PD_Olink_CSF_HC_only %>% filter(str_detect(diagnosis, "HC-ALS"))
PD_HC_list_ALS <- unique(PD_HC_list_ALS$SampleID)


olink_umap_plot(PD_Olink_CSF_HC_only, color_g = "diagnosis", byPanel = FALSE)
olink_pca_plot(PD_Olink_CSF_HC_only, color_g = "diagnosis", byPanel = FALSE)


# PD Clusters
Olink_CSF_PD <- Olink_CSF_combined %>%
        filter(!str_detect(diagnosis, "RBD|HC"))
UMAP_clusters <- olink_umap_plot(Olink_CSF_PD, color_g = "cluster", byPanel = TRUE)

# Outliers
Olink_CSF_outliers <- Olink_CSF_combined %>%
        mutate(outlier = ifelse(SampleID %in% QC_outliers, "Yes", "No"))
UMAP_outliers <- olink_umap_plot(Olink_CSF_outliers, color_g = "outlier", byPanel = TRUE)

# Diagnosis without Outliers 
PD_Olink_CSF_without_outliers_Dx <- PD_Olink_CSF_combined %>%
        mutate(outlier = ifelse(SampleID %in% QC_outliers, "Yes", "No")) %>%
        filter(str_detect(outlier, "No")) %>%
        filter(str_detect(diagnosis, "PD|RBD|HC")) 
UMAP_without_outliers_Dx <- olink_umap_plot(PD_Olink_CSF_without_outliers_Dx, color_g = "diagnosis", byPanel = TRUE)

# Cohorts without Outliers
PD_Olink_CSF_without_outliers_cohorts <- PD_Olink_CSF_combined %>%
        mutate(outlier = ifelse(SampleID %in% QC_outliers, "Yes", "No")) %>%
        filter(str_detect(outlier, "No"))
UMAP_without_outliers_Dx <- olink_umap_plot(PD_Olink_CSF_without_outliers_cohorts, color_g = "cohort", byPanel = TRUE)

# PD clusters without outliers 
Olink_CSF_without_outliers_PD_clusters <- PD_Olink_CSF_combined %>%
        mutate(outlier = ifelse(SampleID %in% QC_outliers, "Yes", "No")) %>%
        filter(str_detect(outlier, "No")) %>%
        filter(!str_detect(diagnosis, "RBD|HC"))
UMAP_without_outliers_PD_clusters <- olink_umap_plot(Olink_CSF_without_outliers_PD_clusters, color_g = "cluster", byPanel = TRUE)

# Plates without outliers
UMAP_without_outliers_plates <- olink_umap_plot(Olink_CSF_without_outliers, color_g = "plate", byPanel = TRUE)

# Gender without outliers 
UMAP_without_outliers_plates <- olink_umap_plot(Olink_CSF_without_outliers, color_g = "gender", byPanel = TRUE)













# import RBD OPDC data and select instruments of interest
OPDC_RBD <- read.csv("~/Documents/R/OPDC/OPDCDiscoveryCohort_ALL_DATA_RBD.csv")

OPDC_RBD_selected <- OPDC_RBD %>% 
        select(subjid, dob, gender, redcap_event_name, x_1_1_rbd_symptom_date, x_1_2_rbd_diag_date, converted_to_pd, conversion_diag_date, conversion_diag_spec, converted_to_pd_r, converted_to_other, converted_to_other_diag,converted_subjid_2, alternate_diag, alternate_diag_1, alternate_diag_specify) %>%
        rename(participant_ID = subjid, date_onset = x_1_1_rbd_symptom_date, date_diagnosis = x_1_2_rbd_diag_date) %>%
        mutate(dob = as.Date(dob)) %>%
        mutate(gender = replace(gender, gender == 0, "male"), gender = replace(gender, gender == 1, "female")) %>%
        rename(visit = redcap_event_name)


# The process of adding dob and gender to all IDs 
OPDC_RBD_selected[OPDC_RBD_selected == ''] <- NA    # changes the character vector " " into NA
OPDC_RBD_dob <- OPDC_RBD_selected[!is.na(OPDC_RBD_selected$dob),]
OPDC_RBD_dob <- subset(OPDC_RBD_dob, select = -c(visit, date_onset, converted_to_pd, converted_to_pd_r, conversion_diag_date, conversion_diag_spec, converted_to_other, converted_to_other_diag,converted_subjid_2, alternate_diag, alternate_diag_1, alternate_diag_specify))

#import list of RBD participants that converted
Michele_list <- read_excel("~/Documents/R/OPDC/OPDC iRBD Converters spreadsheet.xlsx")
Definite_converters <- Michele_list[!is.na(Michele_list$`Confirmed Clinical Date of conversion`),] %>%
        rename(participant_ID = 'ID Number', conversion_dx = 'New Diagnosis', date_of_conversion = 'Confirmed Clinical Date of conversion') %>%
        select(participant_ID, conversion_dx, date_of_conversion)



#Put it all together 
OPDC_RBD_selected_complete <- OPDC_RBD_selected %>%
        subset(select = -c(dob, gender, date_diagnosis, date_onset, converted_to_pd, converted_to_pd_r, conversion_diag_date, conversion_diag_spec, converted_to_other, converted_to_other_diag,converted_subjid_2, alternate_diag, alternate_diag_1, alternate_diag_specify)) %>%
        left_join(OPDC_RBD_dob, by = "participant_ID") %>%
        relocate(dob, gender, date_diagnosis, .after = visit) %>%
        filter(str_detect(visit, "status")) %>%
        left_join(Definite_converters, by = "participant_ID") %>%
        left_join(Olink_CSF_manifest, by = "participant_ID") %>%
        relocate(diagnosis, .after = gender)

#select list of RBD CSF IDs with complete data
RBD_CSF_Olink_tested <- OPDC_RBD_selected_complete %>% 
        filter(str_detect(CSF, "Yes")) %>%
        select(-CSF)







