require(dplyr)
require(stringr)
require(data.table)
require(lubridate)
require(tidyverse)
require(assertthat)
require(readr)

setwd("/pc/resta/RPDR/ProjectForPriya")

rm(list = ls())
data_dir = "data/"
rpdr_file_header = str_c(data_dir, "mns43_100719120614239897_")
rpdr_file_ending = ".txt"
biobank_file_name = str_c(data_dir, "BiobankPortal_mns43_2019-09-18-150014.csv")
output_file_name = str_c(data_dir, "RPDR_cleaned_data.csv")
process_biobankids = TRUE
process_demographics = TRUE
process_diagnoses = TRUE
process_medications = FALSE
process_allergic_rhinitis_meds = FALSE
process_asthma_meds = FALSE
process_labs = FALSE
process_physicalfindings = FALSE
process_deidentified = FALSE

if (process_biobankids){
  BiobankIDs <- data.table(fread(str_c(rpdr_file_header, "Bib", rpdr_file_ending)))
  BiobankIDs <- BiobankIDs %>% rename(Biobank_Subject_ID = Subject_Id) %>% select(Biobank_Subject_ID, EMPI)
  All_merged <- BiobankIDs # now has 2 columns
  rm(BiobankIDs)
}

if (process_demographics){
  Demographics <- data.table(fread(str_c(rpdr_file_header, "Dem", rpdr_file_ending)))
  Demographics <- Demographics %>% select(EMPI, Gender, Date_of_Birth, Age, Date_Of_Death, Race)
  All_merged <- merge(All_merged, Demographics, by = "EMPI") # now has 7 columns (5 new columns)
  rm(Demographics)
}

if (process_diagnoses){
  Diagnoses <- data.table(fread(str_c(rpdr_file_header, "Dia", rpdr_file_ending)))
  
  AAID <- Diagnoses %>% filter(grepl("^[^D].*adren.*insufficiency", Diagnosis_Name)) %>% group_by(EMPI) %>%
    select(EMPI, Diagnosis_Name) %>% unique() %>% summarise(Any_adrenocortical_insufficiency_diagnosis = "Yes")
  All_merged <- left_join(All_merged, AAID, by = "EMPI")
  All_merged <- All_merged %>% mutate(Any_adrenocortical_insufficiency_diagnosis = ifelse(is.na(Any_adrenocortical_insufficiency_diagnosis),
                                                                                          "No",
                                                                                          "Yes"))
  
  CI <- Diagnoses %>% filter(Diagnosis_Name == "Corticoadrenal insufficiency") %>% group_by(EMPI) %>%
    select(EMPI, Diagnosis_Name, Date) %>% unique() %>% mutate(Date = mdy(Date)) %>% arrange(Date) %>%
    summarise(Corticoadrenal_insufficiency_dates = paste(Date, collapse = ";"),
              Corticoadrenal_insufficiency_total_dates = n())
  All_merged <- left_join(All_merged, CI, by = "EMPI")
  
  PAI <- Diagnoses %>% filter(Diagnosis_Name == "Primary adrenocortical insufficiency") %>% group_by(EMPI) %>%
    select(EMPI, Diagnosis_Name, Date) %>% unique() %>% mutate(Date = mdy(Date)) %>% arrange(Date) %>%
    summarise(Primary_adrenocortical_insufficiency_dates = paste(Date, collapse = ";"),
              Primary_adrenocortical_insufficiency_total_dates = n())
  All_merged <- left_join(All_merged, PAI, by = "EMPI")
  
  OAI <- Diagnoses %>% filter(Diagnosis_Name == "Other adrenocortical insufficiency") %>% group_by(EMPI) %>%
    select(EMPI, Diagnosis_Name, Date) %>% unique() %>% mutate(Date = mdy(Date)) %>% arrange(Date) %>%
    summarise(Other_adrenocortical_insufficiency_dates = paste(Date, collapse = ";"),
              Other_adrenocortical_insufficiency_total_dates = n())
  All_merged <- left_join(All_merged, OAI, by = "EMPI")
  
  UAI <- Diagnoses %>% filter(Diagnosis_Name == "Unspecified adrenocortical insufficiency") %>% group_by(EMPI) %>%
    select(EMPI, Diagnosis_Name, Date) %>% unique() %>% mutate(Date = mdy(Date)) %>% arrange(Date) %>%
    summarise(Unspecified_adrenocortical_insufficiency_dates = paste(Date, collapse = ";"),
              Unspecified_adrenocortical_insufficiency_total_dates = n())
  All_merged <- left_join(All_merged, UAI, by = "EMPI")
  
  rm(AAID, CI, OAI, PAI, UAI)
  rm(Diagnoses)
}

if (process_medications){
  Medications <- data.table(fread(str_c(rpdr_file_header, "Med", rpdr_file_ending)))
  
  Bec_search <- str_c("beclomethasone *(-|\\d|i|o|dipropionate($| (-|40.*aer|80.*aer |0|o)))|qvar(\\)|-| \\d.*actuation)",
                    "bec(lovent-|onase i)|vanc(enase i|eril).*oncall", sep = "|")
  Bec <- Medications %>% filter(grepl(Bec_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Beclomethasone_dipropionate_dates = paste(Medication_Date, collapse = ";"),
                                           Beclomethasone_dipropionate_total_dates = n())
  All_merged <- left_join(All_merged, Bec, by = "EMPI")
  
  Bud_search <- "pulmicort.*[^n]$|^budesonide($|/f.*80| (ne|inh(l$|.*powder)|oral s|-|0.*(powder|suspension|ampul|in)|1|9))"
  Bud <- Medications %>% filter(grepl(Bud_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Budesonide_dates = paste(Medication_Date, collapse = ";"),
                                           Budesonide_total_dates = n())
  All_merged <- left_join(All_merged, Bud, by = "EMPI")
  
  Cic_search <- "^ciclesonide|^alvesco 16"
  Cic <- Medications %>% filter(grepl(Cic_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Ciclesonide_dates = paste(Medication_Date, collapse = ";"),
                                           Ciclesonide_total_dates = n())
  All_merged <- left_join(All_merged, Cic, by = "EMPI")
  
  Medications %>% filter(grepl(Cle_search, Medication, ignore.case = TRUE)) %>% group_by(Medication) %>% summarise(med = n())
  
  C_D_Fln <- "ciclesonide (hfa|0)|alvesco 160|dexamethasone sod phospha($|te 4 )|flunisolide(-| (i|\\())|aerobid"
  Flt <- str_c("flovent|(^|id-)fluticasone(-o| (\\d.*(salmetero$|inhaler|blister)|h|inh( |a)",
                    "propionate.*(adap|479(6|7|8)|a)$|furoate \\d))|arnuity|veramyst", sep = "|")
  M_T <- "asmanex|^mometasone (1|2)|azmacort (-|1)|nasacort \\d|^triamcinolone.*(adap|r-oncall)$"
  InCc_search <- str_c(Bec, Bud, C_D_Fln, Flt, M_T, sep = "|")
  InCc <- Medications %>% filter(grepl(InCc_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Inhaled_corticosteroids_dates = paste(Medication_Date, collapse = ";"),
                                           Inhaled_corticosteroids_total_dates = length(Medication_Date))
  All_merged <- left_join(All_merged, InCc, by = "EMPI")
  rm(Bec, Bud, C_D_Fln, Flt, M_T, InCc_search)

  # InC_search <- "nedocromil sodium|^cromolyn (neb|sodium *(-|s)|20)|intal"
  # InC <- Medications %>% filter(grepl(InC_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
  #   select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
  #   arrange(Medication_Date) %>% summarise(Inhaled_chromones_dates = paste(Medication_Date, collapse = ";"),
  #                                          Inhaled_chromones_total_dates = length(Medication_Date))
  # All_merged <- left_join(All_merged, InC, by = "EMPI")
  # rm(InC_search)
  # 
  # InA_U_A_T <- "Umeclidinium.*Blister|Incruse|aclidinium|tudorza|tiotropium|spiriva"
  # InA_I <- str_c("id-ipratropium|^ipratropium (bromide (\\.|0|1).*(solution|m$|aer)|inhaler|nebulizer)|^ipratrop.*(0.02)",
  #                "^ipratropium( |/).*s.*neb|^atrovent(-| (h.*(inhaler|oncall)$|\\d|s))", sep = "|")
  # InA_search <- str_c(InA_U_A_T, InA_I, sep = "|")
  # InA <- Medications %>% filter(grepl(InA_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
  #   select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
  #   arrange(Medication_Date) %>% summarise(Inhaled_anticholinergics_dates = paste(Medication_Date, collapse = ";"),
  #                                          Inhaled_anticholinergics_total_dates = length(Medication_Date))
  # All_merged <- left_join(All_merged, InA, by = "EMPI")
  # rm(InA_U_A_T, InA_I, InA_search)
  # 
  # ICS_LABA_search <- str_c("budesonide.*formotero.*(oncall|actuat)|symbicort|advair|fluticasone.*salmeterol|fluticasone.*vilanterol",
  #                          "breo|mometasone.*formoterol|dulera", sep = "|")
  # ICS_LABA <- Medications %>% filter(grepl(ICS_LABA_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
  #   select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
  #   arrange(Medication_Date) %>% summarise(ICS_LABA_dates = paste(Medication_Date, collapse = ";"),
  #                                          ICS_LABA_total_dates = length(Medication_Date))
  # All_merged <- left_join(All_merged, ICS_LABA, by = "EMPI")
  # rm(ICS_LABA_search)
  # 
  # rm(ICS_LABA, InA, InC, InCc)
  
  rm(Medications)
}

if (process_labs){
  Labs <- data.table(fread(str_c(rpdr_file_header, "Lab", rpdr_file_ending)))
  labs_cats <- Labs %>% group_by(Group_Id) %>% summarise(test = n())
  
  # AEO_S_N <- Labs %>% filter(grepl("AEO", Group_Id), grepl("^(\\d|\\.|<|>)", Result)) %>% group_by(EMPI) %>%
  #   select(EMPI, Group_Id, Result, Seq_Date_Time) %>% unique() %>%
  #   mutate(Seq_Date_Time = mdy_hm(Seq_Date_Time),
  #          Result = gsub("<1", "0.99", Result),
  #          Result = gsub("((\\d|\\.)*).*", "\\1", Result),
  #          Result = as.numeric(Result),
  #          Result = ifelse(Result < 20, Result*1000, Result)) %>%
  #   group_by(EMPI) %>% summarise(AEO_average = mean(Result), AEO_counts = n())
  # All_merged <- left_join(All_merged, AEO_S_N, by = "EMPI")
  # 
  # IGE <- Labs %>% filter(Group_Id == "IGE", grepl("^(\\d|\\.|<|>)", Result)) %>% group_by(EMPI) %>%
  #   select(EMPI, Group_Id, Result, Seq_Date_Time) %>% unique() %>%
  #   mutate(Result = gsub("< *25", "24.9", Result),
  #          Result = gsub("< *15", "10.9", Result),
  #          Result = gsub("< *10", "9.9", Result),
  #          Result = gsub("< *8", "7.9", Result),
  #          Result = gsub("< *7.5", "7.4", Result),
  #          Result = gsub("< *2.*", "1.9", Result),
  #          Result = gsub("< *5", "4.9", Result),
  #          Result = gsub("< *1.5", "1.4", Result),
  #          Result = gsub("> *5000", "5000.1", Result),
  #          Result = as.numeric(Result),
  #          Seq_Date_Time = mdy_hm(Seq_Date_Time)) %>%
  #   group_by(EMPI) %>% summarise(IGE_average = mean(Result), IGE_counts = n())
  # All_merged <- left_join(All_merged, IGE, by = "EMPI") # now has 71 columns
  # 
  # rm(AEO_S_N, IGE)
  rm(Labs)
}

if (process_physicalfindings){
  Phy <- data.table(fread(str_c(rpdr_file_header, "Phy", rpdr_file_ending)))
  phy_cats <- Phy %>% group_by(Concept_Name) %>% summarise(finding = n())
  # FEV1 <- Phy %>% filter(grepl("FEV1", Concept_Name)) %>% unique() %>% group_by(EMPI) %>%
  #   separate(Result, c("Value_Percentage", "Percentage_NA"), sep = " *\\({1,}") %>%
  #   mutate(Percentage_NA = gsub("((\\d|\\.)*).*", "\\1", Percentage_NA),
  #          Value_Percentage = as.numeric(Value_Percentage),
  #          Percentage_NA = as.numeric(Percentage_NA),
  #          Value = ifelse(is.na(Percentage_NA), NA, Value_Percentage),
  #          Percentage = ifelse(is.na(Percentage_NA), Value_Percentage, Percentage_NA)) %>%
  #   summarise(FEV1_value_average = mean(Value, na.rm = TRUE), FEV1_value_counts = sum(!is.na(Value)),
  #             FEV1_percentage_average = mean(Percentage), FEV1_percentage_counts = n()) %>%
  #   mutate(FEV1_value_average = ifelse(is.nan(FEV1_value_average), NA, FEV1_value_average))
  # All_merged <- left_join(All_merged, FEV1, by = "EMPI")
  # 
  # FVC <- Phy %>% filter(grepl("FVC", Concept_Name)) %>% unique() %>% group_by(EMPI) %>%
  #   separate(Result, c("Value_Percentage", "Percentage_NA"), sep = " *\\({1,}") %>%
  #   mutate(Percentage_NA = gsub("((\\d|\\.)*).*", "\\1", Percentage_NA),
  #          Value_Percentage = as.numeric(Value_Percentage),
  #          Percentage_NA = as.numeric(Percentage_NA),
  #          Value = ifelse(is.na(Percentage_NA), NA, Value_Percentage),
  #          Percentage = ifelse(is.na(Percentage_NA), Value_Percentage, Percentage_NA)) %>%
  #   summarise(FVC_value_average = mean(Value, na.rm = TRUE), FVC_value_counts = sum(!is.na(Value)),
  #             FVC_percentage_average = mean(Percentage), FVC_percentage_counts = n()) %>%
  #   mutate(FVC_value_average = ifelse(is.nan(FVC_value_average), NA, FVC_value_average))
  # All_merged <- left_join(All_merged, FVC, by = "EMPI") # now has 79 columns
  # 
  # rm(FEV1, FVC)  
  rm(Phy)
}

if (process_deidentified){
  Deidentified <- read_csv(biobank_file_name)
  names(Deidentified) <- gsub("(\\(|\\))", "",
                              gsub("_-_", "_",
                                   gsub("2._(.*)", "\\1",
                                        gsub(" ", "_",
                                             gsub("(.*)( \\[.*\\])", "\\1", names(Deidentified))))))
  Deidentified <- Deidentified %>% select(-c(Patient_Number, People_with_genomic_data)) %>%
    rename(Asthma = Asthma_current_or_past_history_PPV_0.90,
           Control = `0-_10-year_survival_probability_=_98.30%`,
           Highest_education = `What_is_the_highest_grade_in_school_that_you_finished?`,
           Smoked_100_cigarettes = `Have_you_smoked_at_least_100_cigarettes_in_your_lifetime?`) %>%
    mutate(Asthma = ifelse(Asthma == "Yes", 1, 0),
           Control = ifelse(Control == "Yes", 1, 0),
           BMI = as.numeric(BMI),
           Age_at_Survey = as.numeric(gsub("\\[(.*)\\]", "\\1", Age_at_Survey)),
           Highest_education = ifelse(is.na(Highest_education), NA,
                                      gsub("\\[ \\d. (.*)\\]", "\\1", Highest_education)),
           Smoked_100_cigarettes = ifelse(is.na(Smoked_100_cigarettes), NA,
                                          ifelse(grepl("^\\[ 3. No,( |\\w)*\\]$", Smoked_100_cigarettes), 0, 1)))

  All_merged <- left_join(All_merged, Deidentified, by = "Biobank_Subject_ID") # now has 86
  All_merged <- All_merged %>% mutate(Allergic_rhinitis = ifelse(is.na(Allergic_rhinitis_first_diagnosis), 0, 1)) %>%
    filter(!(Allergic_rhinitis == 1 & Control == 1),
           !(Asthma == 1 & Control == 1),
           !(Allergic_rhinitis == 0 & Asthma == 0 & Control == 0)) %>%
    select(EMPI, Biobank_Subject_ID, Allergic_rhinitis, Asthma, Control,
           Gender, Date_of_Birth, Date_Of_Death, Age, Age_at_Survey, Race,
           Smoked_100_cigarettes, Highest_education, BMI, AEO_average:FVC_percentage_counts,
           Allergic_rhinitis_first_diagnosis:Xanthines_total_dates)
  
  rm(Deidentified)
}

fwrite(All_merged, output_file_name)
