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
process_medications = TRUE
process_labs = FALSE
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
  All_merged <- All_merged %>% mutate(Any_adrenocortical_insufficiency_diagnosis =
                                        ifelse(is.na(Any_adrenocortical_insufficiency_diagnosis),
                                               "No", "Yes"))
  
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
  Bud_search <- "pulmicort.*[^n]$|^budesonide($|/f.*80| (ne|inh(l$|.*powder)|oral s|-|0.*(powder|suspension|ampul|in)|1|9))"
  Cic_search <- "^Ciclesonide|^Alvesco 16"
  Dex_search <- "Dexamethasone sod phospha($|te 4 )"
  Fln_search <- "Flunisolide(-| (i|\\())|Aerobid"
  Flt_search <- str_c("Flovent|(^|id-)Fluticasone(-o| (\\d.*(salmetero$|inhaler|blister)|h|inh( |a)",
                      "propionate.*(adap|479(6|7|8)|a)$|furoate \\d))|arnuity|veramyst", sep = "|")
  Flt_sal_search <- "advair|fluticasone.*salmeterol"
  Mom_search <- "Asmanex|^Mometasone (1|2)"
  Tri_search <- "Azmacort (-|1)|Nasacort \\d|^Triamcinolone.*(adap|r-oncall)$"
  Any_search <- str_c(Bec_search, Bud_search, Cic_search, Dex_search, Fln_search, Flt_search,
                      Flt_sal_search, Mom_search, Tri_search, sep = "|")
  
  AAB <- Medications %>% filter(grepl(Any_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication) %>% unique() %>% summarise(Any_antiasthma_bronchodialator_prescription = "Yes")
  All_merged <- left_join(All_merged, AAB, by = "EMPI")
  All_merged <- All_merged %>% mutate(Any_antiasthma_bronchodialator_prescription =
                                        ifelse(is.na(Any_antiasthma_bronchodialator_prescription),
                                               "No", "Yes"))
  
  Bec <- Medications %>% filter(grepl(Bec_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Beclomethasone_dipropionate_dates = paste(Medication_Date, collapse = ";"),
                                           Beclomethasone_dipropionate_total_dates = n())
  All_merged <- left_join(All_merged, Bec, by = "EMPI")
  
  Bud <- Medications %>% filter(grepl(Bud_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Budesonide_dates = paste(Medication_Date, collapse = ";"),
                                           Budesonide_total_dates = n())
  All_merged <- left_join(All_merged, Bud, by = "EMPI")
  
  Cic <- Medications %>% filter(grepl(Cic_search, Medication)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Ciclesonide_dates = paste(Medication_Date, collapse = ";"),
                                           Ciclesonide_total_dates = n())
  All_merged <- left_join(All_merged, Cic, by = "EMPI")
  
  Dex <- Medications %>% filter(grepl(Dex_search, Medication)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Dexamethasone_dates = paste(Medication_Date, collapse = ";"),
                                           Dexamethasone_total_dates = n())
  All_merged <- left_join(All_merged, Dex, by = "EMPI")
  
  Fln <- Medications %>% filter(grepl(Fln_search, Medication)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Flunisolide_dates = paste(Medication_Date, collapse = ";"),
                                           Flunisolide_total_dates = n())
  All_merged <- left_join(All_merged, Fln, by = "EMPI")
  
  Flt <- Medications %>% filter(grepl(Flt_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Fluticasone_dates = paste(Medication_Date, collapse = ";"),
                                           Fluticasone_total_dates = n())
  All_merged <- left_join(All_merged, Flt, by = "EMPI")
  
  Flt_sal <- Medications %>% filter(grepl(Flt_sal_search, Medication, ignore.case = TRUE)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Fluticasone_salmeterol_dates = paste(Medication_Date, collapse = ";"),
                                           Fluticasone_salmeterol_total_dates = n())
  All_merged <- left_join(All_merged, Flt_sal, by = "EMPI")

  Mom <- Medications %>% filter(grepl(Mom_search, Medication)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Mometasone_dates = paste(Medication_Date, collapse = ";"),
                                           Mometasone_total_dates = n())
  All_merged <- left_join(All_merged, Mom, by = "EMPI")
  
  Tri <- Medications %>% filter(grepl(Tri_search, Medication)) %>% group_by(EMPI) %>%
    select(EMPI, Medication_Date) %>% unique() %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(Medication_Date) %>% summarise(Triamcinolone_dates = paste(Medication_Date, collapse = ";"),
                                           Triamcinolone_total_dates = n())
  All_merged <- left_join(All_merged, Tri, by = "EMPI")
  
  # Medications %>% filter(grepl(Tri_search, Medication, ignore.case = TRUE)) %>% group_by(Medication) %>% summarise(med = n())
  rm(Any_search, Bec_search, Bud_search, Cic_search, Dex_search, Fln_search, Flt_search,
     Flt_sal_search, Mom_search, Tri_search)
  rm(AAB, Bec, Bud, Cic, Dex, Fln, Flt, Flt_sal, Mom, Tri)
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
