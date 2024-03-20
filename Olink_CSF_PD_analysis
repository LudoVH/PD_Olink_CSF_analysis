

if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "lubridate",
               "OlinkAnalyze",
               "rio")


###################### EXPLORATION OF FINAL DATA SET ########################

# Import QC'ed data (ALS + PD cohorts)
PD_Olink_CSF_data_clean_without_bad_samples <- import("PD_Olink_CSF_data_clean_without_bad_samples.csv")

# PD cohort data only         
PD_Olink_CSF_data_PD_cohort <- PD_Olink_CSF_data_clean_without_bad_samples %>%
        filter(cohort == "PD") %>%
        select(-QC_WARN_percentage, -QC_rate) %>%
        mutate(date_onset = year(date_onset), date_diagnosis = year(date_diagnosis), csf_date = year(csf_date))



export(PD_Olink_CSF_data_PD_cohort, "PD_Olink_CSF_data_PD_cohort.csv")

# UMAP labeling 
PD_Olink_CSF_analysis_disease_duration <- PD_Olink_CSF_data_PD_cohort_single_visit %>%
        d


# PD + RBD vs HC 
PD_Olink_CSF_PD_RBD_vs_HC <- PD_Olink_CSF_data_PD_cohort %>%
        mutate(asyn = ifelse(diagnosis == "RBD" | diagnosis == "PD", "a-Syn", "control"))

PD_Olink_CSF_PD_RBD_vs_HC_table <- PD_Olink_CSF_PD_RBD_vs_HC %>%
        distinct(SampleID, .keep_all = TRUE) %>%
        tbl_summary(include = diagnosis)
print(PD_Olink_CSF_PD_RBD_vs_HC_table)

# perform t-test
ttest_results_PD_Olink_CSF_PD_RBD_vs_HC <- olink_ttest(df = PD_Olink_CSF_PD_RBD_vs_HC,
                                                       variable = 'asyn')
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_PD_Olink_CSF_PD_RBD_vs_HC,
                   x_lab = 'a-Syn disease', olinkid_list = FALSE)


# PD vs HC
PD_Olink_CSF_PD_vs_HC <- PD_Olink_CSF_data_PD_cohort %>%
        filter(diagnosis == "PD" | diagnosis == "HC")

# perform t-test
ttest_results_PD_vs_HC <- olink_ttest(df = PD_Olink_CSF_PD_vs_HC,
                                      variable = 'diagnosis')

# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_PD_vs_HC,
                   x_lab = 'Disease status', olinkid_list = FALSE)


# RBD v HC
PD_Olink_CSF_RBDvHC <- PD_Olink_CSF_data_PD_cohort %>%
        filter(str_detect(diagnosis,"RBD|HC"))

# perform t-test
ttest_results_RBDvHC <- olink_ttest(df = PD_Olink_CSF_RBDvHC,
                                    variable = 'diagnosis')

# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results_RBDvHC,
                   x_lab = 'RBD vs HC')


# PD Clusters
PD_Olink_CSF_clusters <- PD_Olink_CSF_data_final %>%
        filter(!str_detect(diagnosis, "RBD|HC")) %>%
        mutate(cluster = ifelse(is.na(cluster), "No cluster", cluster)) %>%
        filter(!cluster == "No cluster")
olink_umap_plot(PD_Olink_CSF_clusters, color_g = "cluster", byPanel = FALSE)

PCA_clusters <- olink_pca_plot(PD_Olink_CSF_clusters, color_g = "cluster", byPanel = TRUE)













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
