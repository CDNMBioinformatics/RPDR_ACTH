require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(data.table) # fread, fwrite
require(lubridate) # date time conversion
require(tidyverse)
require(logging)
require(yaml)

################################
###### Create main dataset #####
################################

rm(list = ls())
timestamp <- format(Sys.time(), "%Y-%m-%d")
config <- yaml.load_file("configDosageQuery.yml")
output_file_name = str_c(config$output_file_header, timestamp, config$general_file_ending)
if (config$create_intermediates){
  if (!dir.exists(config$intermediate_files_dir)){
    dir.create(config$intermediate_files_dir, recursive = TRUE)
  }
}
if (!dir.exists(config$log_params$path)){
  dir.create(config$log_params$path, recursive = TRUE)
}
config$log_params$file_name <- str_c(config$log_params$path, timestamp, config$log_params$file_ending)
logReset()
setLevel(config$log_params$level, getLogger())
addHandler(writeToConsole, level = config$log_params$level)
addHandler(writeToFile, file = config$log_params$file_name, level = config$log_params$level)

source(str_c(config$rpdr_src_dir, "Start_Cleaning.R")) # start_processing
source(str_c(config$rpdr_src_dir, "Add_Biobank_File.R")) # process_deidentified
source(str_c(config$rpdr_src_dir, "Add_Plasma_File.R")) # process_plasma
source(str_c(config$rpdr_src_dir, "Clean_Health_History_File.R")) # process_physical_set_range
source(str_c(config$rpdr_src_dir, "Clean_Diagnoses_File.R")) # process_diagnoses_set_range
source(str_c(config$rpdr_src_dir, "Clean_Medications_File.R"))
source(str_c(config$rpdr_src_dir, "Clean_Labs_File.R"))


# Begin processing and cleaning of RPDR files: 10 columns should be created
All_merged <- start_processing()

# Should add 5 columns -> 15 columns
All_merged <- process_deidentified(clean.list = TRUE)

# Should add 4 columns -> 19 columns
All_merged <- process_plasma(date_cutoff = "2020-06-30", reduce_to_plasma_only = TRUE)
# Add 1 more column so getting the range of data is easier-> 20 columns
All_merged <- All_merged %>% mutate(Plasma_First_Minus_5_Years = First_Collection_Date %m-% years(5))

# Get file of Genotype Ids 2 more columns -> 22
Genotype_File <- fread(config$genotype_file_name)
Genotype_File <- Genotype_File %>% rename("Biobank_Subject_ID" = SubjectID)
All_merged <- left_join(All_merged, Genotype_File) %>%
  mutate(Has_Genotype_Info_MegE = ifelse(is.na(SentrixID), "0", "1")) %>%
  select(EMPI:All_Collection_Dates, Has_Genotype_Info_MegE, SentrixID, everything())
rm(Genotype_File)

# Should add 28 columns -> 50 columns
All_merged <- process_physical_set_range(Return_Influenza = FALSE,
                                         min_dates = "Plasma_First_Minus_5_Years",
                                         max_dates = "First_Collection_Date")

# Should add 26 (6 + 4*5) columns -> 76 columns
All_merged <- process_diagnoses_set_range(Diagnoses_Of_Interest =
                                            list("Adrenal insufficiency" =
                                                   c("Corticoadrenal insufficiency",
                                                     "Primary adrenocortical insufficiency",
                                                     "Other adrenocortical insufficiency",
                                                     "Unspecified adrenocortical insufficiency")),
                                          min_dates = "Plasma_First_Minus_5_Years",
                                          max_dates = "First_Collection_Date")

# Should add 36 (5*6 + 6) columns -> 112 columns
All_merged <- process_diagnoses_set_range(Diagnoses_Of_Interest =
                                            list("Mild intermittent asthma" = c("Extrinsic asthma",
                                                                                "Intrinsic asthma",
                                                                                "Mild intermittent asthma"),
                                                 "Mild persistent asthma" = c("Mild persistent asthma"),
                                                 "Moderate persistent asthma" = c("Moderate persistent"),
                                                 "Severe persistent asthma" = c("Severe persistent asthma"),
                                                 "Other and unspecified asthma" = c("^Asthma(, acute|-|$)",
                                                                                    "Asthma,{0,1} unspecified",
                                                                                    "Asthmatic bronchitis",
                                                                                    "Cough variant asthma",
                                                                                    "Exercise.induced asthma",
                                                                                    "Other asthma",
                                                                                    "Unspecified asthma")),
                                          Exact = FALSE,
                                          Individual_Info = FALSE,
                                          Merge_Group_Info_Name = "Asthma",
                                          min_dates = "Plasma_First_Minus_5_Years",
                                          max_dates = "First_Collection_Date")

# Should add 12 (2*6) columns -> 124 columns
All_merged <- process_diagnoses_set_range(Diagnoses_Of_Interest =
                                            list("Bronchiestasis" = c("^Bronchiectasis"),
                                                 "Other chronic obstructive pulmonary disease" =
                                                   c("Chronic airway",
                                                     "Chronic bronchitis-",
                                                     "Chronic obstructive asthma",
                                                     "Chronic obstructive lung",
                                                     "Chronic obstructive pulmonary disease.[^\\(]",
                                                     "Obstructive chronic bronchitis")),
                                          Exact = FALSE,
                                          Individual_Info = FALSE,
                                          min_dates = "Plasma_First_Minus_5_Years",
                                          max_dates = "First_Collection_Date")

# Should add 48 (7*6 + 6) columns -> 172 columns
All_merged <- process_diagnoses_set_range(Diagnoses_Of_Interest =
                                            list("Mild intermittent asthma with acute exacerbation" =
                                                   c("Extrinsic asthma with acute exacerbation",
                                                     "Intrinsic asthma, with acute exacerbation",
                                                     "Mild intermittent asthma with \\(acute\\) exacerbation"),
                                                 "Mild persistent asthma with acute exacerbation" =
                                                   c("Mild persistent asthma with \\(acute\\) exacerbation"),
                                                 "Moderate persistent asthma with acute exacerbation" =
                                                   c("Moderate persistent asthma with \\(acute\\) exacerbation"),
                                                 "Severe persistent asthma with acute exacerbation" =
                                                   c("Severe persistent asthma with \\(acute\\) exacerbation"),
                                                 "Other and unspecified asthma with acute exacerbation" =
                                                   c("Asthma, acute exacerbation-LMR 1288",
                                                     "Asthma, unspecified type, with acute exacerbation",
                                                     "Unspecified asthma with \\(acute\\) exacerbation"),
                                                 "Bronchiectasis with acute exacerbation" =
                                                   c("Bronchiectasis with \\(acute\\) exacerbation",
                                                     "Bronchiectasis with acute exacerbation"),
                                                 "Other chronic obstructive pulmonary disease with acute exacerbation" =
                                                   c("Chronic obstructive asthma with acute exacerbation",
                                                     "Chronic obstructive pulmonary disease with \\(acute\\) exacerbation",
                                                     "Obstructive chronic bronchitis with acute exacerbation")),
                                          Exact = FALSE,
                                          Individual_Info = FALSE,
                                          Merge_Group_Info_Name = "Acute exacerbation",
                                          min_dates = "Plasma_First_Minus_5_Years",
                                          max_dates = "First_Collection_Date")

# Should add 373 (9 + 13*28) columns -> 545 columns
All_merged <- process_medications_set_range(Medications_Of_Interest =
                                              list("Anti-inflammatories, inhalation" =
                                                     c("Beclomethasone dipropionate",
                                                       "Budesonide",
                                                       "Ciclesonide",
                                                       "Dexamethasone",
                                                       "Flunisolide",
                                                       "Fluticasone",
                                                       "Mometasone",
                                                       "Triamcinolone"),
                                                   "Antiasthma, other" =
                                                     c("Budesonide/formoterol",
                                                       "Fluticasone/salmeterol",
                                                       "Fluticasone/vilanterol",
                                                       "Formoterol/mometasone"),
                                                   "Glucocorticoids" =
                                                     c("Betamethasone")),
                                            Group_Info = FALSE,
                                            merged_group_name = "Inhaled Corticosteroids",
                                            min_dates = "Plasma_First_Minus_5_Years",
                                            max_dates = "First_Collection_Date",
                                            Daily_Dose_Info = TRUE)

# Should add 142 (9 + 19*7) columns -> 687 columns
All_merged <- process_medications_set_range(Medications_Of_Interest =
                                              list("Antiasthma, antileukotrienes" =
                                                     c("Ipratropium/albuterol"),
                                                   "Antiasthma, other" =
                                                     c("Albuterol/ipratropium",
                                                       "Montelukast",
                                                       "Zafirlukast",
                                                       "Zileuton"),
                                                   "Bronchodilators, anticholinergic" =
                                                     c("Ipratropium",
                                                       "Tiotropium"),
                                                   "Bronchodilators, sympathomimetic, inhalation" =
                                                     c("Albuterol",
                                                       "Levalbuterol"),
                                                   "Bronchodilators, sympathomimetic, oral" =
                                                     c("Albuterol"),
                                                   "Glucocorticoids" =
                                                     c("Budesonide",
                                                       "Cortisone",
                                                       "Dexamethasone",
                                                       "Hydrocortisone",
                                                       "Methylprednisolone",
                                                       "Prednisolone",
                                                       "Prednisone",
                                                       "Triamcinolone"),
                                                   "Immune suppressants" =
                                                     c("Omalizumab")),
                                            Group_Info = FALSE,
                                            merged_group_name = "Not Inhaled Corticosteroids",
                                            min_dates = "Plasma_First_Minus_5_Years",
                                            max_dates = "First_Collection_Date")

# Should add 189 columns -> 876 columns
All_merged <- process_ACTH_labs_set_range(min_dates = "Plasma_First_Minus_5_Years",
                                          max_dates = "First_Collection_Date")
# Should add 189 (9*21*21) columns -> 1065 columns
All_merged <- process_ACTH_labs_set_range(strict = TRUE,
                                          min_dates = "Plasma_First_Minus_5_Years",
                                          max_dates = "First_Collection_Date")

fwrite(All_merged, output_file_name)
loginfo("Processing Complete")

# 
# ########################################
# ###### Create additional data sets #####
# ########################################
# config$bonus_dir <- str_c(config$data_dir, "AdditionalDataSets/")
# if (!dir.exists(config$bonus_dir)){ dir.create(config$bonus_dir, recursive = TRUE) }
# 
# # timestamp <- format(Sys.time() - 86400, "%Y-%m-%d")
# # input_file_name = str_c(config$output_file_header, timestamp, config$general_file_ending)
# Cleaned_name <- str_c("RPDR_cleaned_data")
# # RPDR_cleaned_data <- fread(input_file_name)
# assign(Cleaned_name, All_merged)
# 
# Summary_DF <- data.frame(File_Name = Cleaned_name,
#                          nSubjects = nrow(get(Cleaned_name)),
#                          Description = "Cleaned asthma dataset")
# DF_names <- c(Cleaned_name)
# 
# #(1) PPV 0.85
# {
#   PPV_0.85 <- get(Cleaned_name) %>% filter(Asthma_current_or_past_history_custom_PPV_0_80PPV_List_of_All_Values_To_6_30_2020 >= 0.85)
#   loginfo(str_c(nrow(PPV_0.85), " subjects available with PPV >= 0.85"))
#   fwrite(PPV_0.85, str_c(config$bonus_dir, "PPV_0.85_", timestamp, config$general_file_ending))
#   Summary_DF <- rbind(Summary_DF,
#                       data.frame(File_Name = "PPV_0.85",
#                                  nSubjects = nrow(PPV_0.85),
#                                  Description = "Subjects available with PPV >= 0.85"))
#   DF_names <- c(DF_names, deparse(substitute(PPV_0.85)))
# }
# 
# #(2) ICS use
# ICS_cutoffs <- c(1, 4, 10)
# for (ic in ICS_cutoffs){
#   ic_name = ifelse(ic == 1, "ICS", str_c("ICS_", ic))
#   assign(ic_name, PPV_0.85 %>% filter(Any_Inhaled_Corticosteroids_total_dates >= ic))
#   loginfo(str_c(nrow(get(ic_name)), " subjects have at least ", ic," ICS prescription(s)"))
#   fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
#   Summary_DF <- rbind(Summary_DF,
#                       data.frame(File_Name = ic_name,
#                                  nSubjects = nrow(get(ic_name)),
#                                  Description = str_c("Subjects available with PPV >= 0.85",
#                                                      " and have at least ", ic," ICS prescription(s)")))
#   DF_names <- c(DF_names, ic_name)
# }
# rm(ic, ic_name)
# 
# #(*) Exacerbations
# for (ic in ICS_cutoffs){
#   ic_name = ifelse(ic == 1, "ICS", str_c("ICS_", ic))
#   ic_c_name = str_c(ic_name, "_Exacerbation")
#   assign(ic_c_name, get(ic_name) %>% filter(Any_Acute_exacerbation == "Yes"))
#   loginfo(str_c(nrow(get(ic_c_name)), " subjects have at least ", ic," ICS prescription(s) and ",
#                 "had at least 1 exacerbation diagnosis"))
#   fwrite(get(ic_c_name), str_c(config$bonus_dir, ic_c_name, "_", timestamp, config$general_file_ending))
#   Summary_DF <- rbind(Summary_DF,
#                       data.frame(File_Name = ic_c_name,
#                                  nSubjects = nrow(get(ic_c_name)),
#                                  Description = str_c("Subjects available with PPV >= 0.85",
#                                                      " have at least ", ic," ICS prescription(s)",
#                                                      " and have at least 1 exacerbation diagnosis")))
#   DF_names <- c(DF_names, ic_c_name)
# }
# rm(ic, ic_name, ic_c_name)
# 
# exacerbation_diagnoses <- c(5, 10, 15, 20)
# for (cc in exacerbation_diagnoses){
#   for (ic in ICS_cutoffs){
#     ic_c_name = str_c(ifelse(ic == 1, "ICS", str_c("ICS_", ic)), "_Exacerbation")
#     ic_c_c_name = str_c(ic_c_name, "_Count_", cc)
#     assign(ic_c_c_name, get(ic_c_name) %>% filter(Any_Acute_exacerbation_total_diagnoses >= cc))
#     loginfo(str_c(nrow(get(ic_c_c_name)), " subjects have at least ", ic," ICS prescription(s), ",
#                   "and at least ", cc, " exacerbation diagnoses"))
#     fwrite(get(ic_c_c_name), str_c(config$bonus_dir, ic_c_c_name, "_", timestamp, config$general_file_ending))
#     Summary_DF <- rbind(Summary_DF,
#                         data.frame(File_Name = ic_c_c_name,
#                                    nSubjects = nrow(get(ic_c_c_name)),
#                                    Description = str_c("Subjects available with PPV >= 0.85",
#                                                        " have at least ", ic," ICS prescription(s)",
#                                                        " and have at least ", cc, " exacerbation diagnoses")))
#     DF_names <- c(DF_names, ic_c_c_name)
#   }
# }
# rm(cc, ic, ic_c_name, ic_c_c_name)
# 
# 
# 
# #(3a) ICS use and cortisol test
# for (ic in ICS_cutoffs){
#   ic_name = ifelse(ic == 1, "ICS", str_c("ICS_", ic))
#   ic_c_name = str_c(ic_name, "_Tested")
#   assign(ic_c_name, get(ic_name) %>% filter(All_Cortisol_nTotalDates >= 1))
#   loginfo(str_c(nrow(get(ic_c_name)), " subjects have at least ", ic," ICS prescription(s) and ",
#                 "had at least 1 cortisol test"))
#   fwrite(get(ic_c_name), str_c(config$bonus_dir, ic_c_name, "_", timestamp, config$general_file_ending))
#   Summary_DF <- rbind(Summary_DF,
#                       data.frame(File_Name = ic_c_name,
#                                  nSubjects = nrow(get(ic_c_name)),
#                                  Description = str_c("Subjects available with PPV >= 0.85",
#                                                      " have at least ", ic," ICS prescription(s)",
#                                                      " and have at least 1 cortisol test")))
#   DF_names <- c(DF_names, ic_c_name)
# }
# rm(ic, ic_name, ic_c_name)
# 
# #(3b) ICS use and cortisol test with cutoffs
# cortisol_cutoffs <- c(2, 5, 10, 15, 20)
# for (cc in cortisol_cutoffs){
#   for (ic in ICS_cutoffs){
#     ic_c_name = str_c(ifelse(ic == 1, "ICS", str_c("ICS_", ic)), "_Tested")
#     ic_c_c_name = str_c(ic_c_name, "_Cutoff_", cc)
#     assign(ic_c_c_name, get(ic_c_name) %>% filter(All_Cortisol_Overall_Median_Result <= cc))
#     loginfo(str_c(nrow(get(ic_c_c_name)), " subjects have at least ", ic," ICS prescription(s), ",
#                   "had at least 1 cortisol test, and the median value was no greater than ", cc))
#     fwrite(get(ic_c_c_name), str_c(config$bonus_dir, ic_c_c_name, "_", timestamp, config$general_file_ending))
#     Summary_DF <- rbind(Summary_DF,
#                         data.frame(File_Name = ic_c_c_name,
#                                    nSubjects = nrow(get(ic_c_c_name)),
#                                    Description = str_c("Subjects available with PPV >= 0.85",
#                                                        " have at least ", ic," ICS prescription(s)",
#                                                        " and have at least 1 cortisol test",
#                                                        " where the median value is not greater than ", cc)))
#     DF_names <- c(DF_names, ic_c_c_name)
#   }
# }
# rm(cc, ic, ic_c_name, ic_c_c_name)
# 
# 
# #(4) ICS use and cortisol test and adrenal diagnosis
# for (ic in ICS_cutoffs){
#   ic_c_name = str_c(ifelse(ic == 1, "ICS", str_c("ICS_", ic)), "_Tested")
#   ic_c_d_name = str_c(ic_c_name, "_Diagnosis")
#   assign(ic_c_d_name, get(ic_c_name) %>% filter(Any_Adrenal_insufficiency == "Yes"))
#   loginfo(str_c(nrow(get(ic_c_d_name)), " subjects have at least ", ic," ICS prescription(s) and ",
#                 "had at least 1 cortisol test, and at least 1 adrenal insufficiency diagnosis"))
#   fwrite(get(ic_c_d_name), str_c(config$bonus_dir, ic_c_d_name, "_", timestamp, config$general_file_ending))
#   Summary_DF <- rbind(Summary_DF,
#                       data.frame(File_Name = ic_c_d_name,
#                                  nSubjects = nrow(get(ic_c_d_name)),
#                                  Description = str_c("Subjects available with PPV >= 0.85",
#                                                      " have at least ", ic," ICS prescription(s)",
#                                                      " have at least 1 cortisol test,",
#                                                      " and have at least 1 adrenal insufficency diagnosis")))
#   DF_names <- c(DF_names, ic_c_d_name)
# }
# rm(ic, ic_c_name, ic_c_d_name)
# 
# fwrite(Summary_DF, str_c(config$bonus_dir, "Summary_Overview", config$general_file_ending))
# 
# 
# #(5) ICS changes before and after cortisol test
# # (5a) Compare Medications
# # Should add 183 (2*(9 + 9*7) + 39) columns -> 841 columns
# ICS_Tested_Compare <- process_medications_date_compare_cutoff(Medications_Of_Interest =
#                                                list("Anti-inflammatories, inhalation" =
#                                                       c("Beclomethasone dipropionate", "Budesonide",
#                                                         "Ciclesonide", "Dexamethasone",
#                                                         "Flunisolide", "Fluticasone",
#                                                         "Mometasone", "Triamcinolone"),
#                                                     "Antiasthma, other" = c("Fluticasone/salmeterol")),
#                                              Group_Info = FALSE,
#                                              merged_group_name = "Inhaled Corticosteroids",
#                                              cutoff_variable = deparse(substitute(All_Cortisol_Overall_Median_Result_Date_First_or_closest_below)))
# 
# ic_name = "ICS_preTest"
# assign(ic_name, ICS_Tested_Compare %>% filter(Compare_Any_Inhaled_Corticosteroids == "Yes - Predates cutoff"))
# loginfo(str_c(nrow(get(ic_name)), " subjects were only prescribed any ICS prescription BEFORE their first cortisol test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " and they were only prescribed ICS medications before their test")))
# DF_names <- c(DF_names, ic_name)
# 
# ic_name = "ICS_postTest"
# assign(ic_name, ICS_Tested_Compare %>% filter(Compare_Any_Inhaled_Corticosteroids == "Yes - Postdates cutoff"))
# loginfo(str_c(nrow(get(ic_name)), " subjects were only prescribed any ICS prescription AFTER their first cortisol test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " and they were only prescribed ICS medications after their test")))
# DF_names <- c(DF_names, ic_name)
# 
# ic_name = "ICS_prepostTest_nochange"
# assign(ic_name, ICS_Tested_Compare %>% filter(Compare_Any_Inhaled_Corticosteroids == "Yes - Both ranges",
#                                               Compare_Any_Inhaled_Corticosteroids_most_common_prescription_type == "Same"))
# loginfo(str_c(nrow(get(ic_name)), " subjects had no change in their most common type of ICS prescription before and after their first cortisol test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " and there was no change in the most common ICS prescription type",
#                                                    " before  and after their first cortisol test")))
# DF_names <- c(DF_names, ic_name)
# 
# ICS_ChangedTest <- ICS_Tested_Compare %>% filter(Compare_Any_Inhaled_Corticosteroids_most_common_prescription_type == "Different")
# ic_name = "ICS_prepostTest_difference"
# assign(ic_name, ICS_Tested_Compare %>% filter(Compare_Any_Inhaled_Corticosteroids == "Yes - Both ranges",
#                                               Compare_Any_Inhaled_Corticosteroids_most_common_prescription_type == "Different"))
# loginfo(str_c(nrow(get(ic_name)), " subjects changed their most common type of ICS prescription before and after their first cortisol test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " and changed the most common ICS prescription type",
#                                                    " before and after their first cortisol test")))
# DF_names <- c(DF_names, ic_name)
# 
# # (5b) Compare Acute Exacerbation Diagnoses
# ICS_Tested_Compare <- process_diagnoses_date_compare_cutoff(DF_to_fill = ICS_Tested_Compare,
#                                            Diagnoses_Of_Interest =
#                                              list("Mild intermittent asthma with acute exacerbation" =
#                                                     c("Extrinsic asthma with acute exacerbation",
#                                                       "Intrinsic asthma, with acute exacerbation",
#                                                       "Mild intermittent asthma with \\(acute\\) exacerbation"),
#                                                   "Mild persistent asthma with acute exacerbation" =
#                                                     c("Mild persistent asthma with \\(acute\\) exacerbation"),
#                                                   "Moderate persistent asthma with acute exacerbation" =
#                                                     c("Moderate persistent asthma with \\(acute\\) exacerbation"),
#                                                   "Severe persistent asthma with acute exacerbation" =
#                                                     c("Severe persistent asthma with \\(acute\\) exacerbation"),
#                                                   "Other and unspecified asthma with acute exacerbation" = 
#                                                     c("Asthma, acute exacerbation-LMR 1288",
#                                                       "Asthma, unspecified type, with acute exacerbation",
#                                                       "Unspecified asthma with \\(acute\\) exacerbation"),
#                                                   "Bronchiectasis with acute exacerbation" =
#                                                     c("Bronchiectasis with \\(acute\\) exacerbation",
#                                                       "Bronchiectasis with acute exacerbation"),
#                                                   "Other chronic obstructive pulmonary disease with acute exacerbation" =
#                                                     c("Chronic obstructive asthma with acute exacerbation",
#                                                       "Chronic obstructive pulmonary disease with \\(acute\\) exacerbation",
#                                                       "Obstructive chronic bronchitis with acute exacerbation")),
#                                            Exact = FALSE,
#                                            Individual_Info = FALSE,
#                                            Merge_Group_Info_Name = "Acute exacerbation",
#                                            cutoff_variable = deparse(substitute(All_Cortisol_Overall_Median_Result_Date_First_or_closest_below)))
# 
# ic_name = "Exacerbation"
# assign(ic_name, ICS_Tested_Compare %>% filter(Any_Acute_exacerbation == "Yes"))
# loginfo(str_c(nrow(get(ic_name)), " subjects had at least 1 ICS prescription, at least 1 cortisol test, and at least 1 acute exacerbation diagnosis"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " and at least 1 acute exacerbation diagnosis")))
# DF_names <- c(DF_names, ic_name)
# 
# ic_name = "Exacerbation_preTest"
# assign(ic_name, ICS_Tested_Compare %>% filter(Any_Acute_exacerbation == "Yes",
#                                               Compare_Any_Acute_exacerbation == "Yes - Predates cutoff"))
# loginfo(str_c(nrow(get(ic_name)), " subjects had at least 1 ICS prescription, at least 1 cortisol test, but were only diagnosed with any acute exacerbation before their first test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " but were only diagnosed with any acute exacerbation before their first cortisol test")))
# DF_names <- c(DF_names, ic_name)
# 
# ic_name = "Exacerbation_postTest"
# assign(ic_name, ICS_Tested_Compare %>% filter(Any_Acute_exacerbation == "Yes",
#                                               Compare_Any_Acute_exacerbation == "Yes - Postdates cutoff"))
# loginfo(str_c(nrow(get(ic_name)), " subjects had at least 1 ICS prescription, at least 1 cortisol test, but were diagnosed with any acute exacerbation before and after their first cortisol test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " but were only diagnosed with any acute exacerbation after thir first cortisol tests")))
# DF_names <- c(DF_names, ic_name)
# 
# ic_name = "Exacerbation_prepostTest"
# assign(ic_name, ICS_Tested_Compare %>% filter(Any_Acute_exacerbation == "Yes",
#                                               Compare_Any_Acute_exacerbation == "Yes - Both ranges"))
# loginfo(str_c(nrow(get(ic_name)), " subjects had at least 1 ICS prescription, at least 1 cortisol test, and were diagnosed with any acute exacerbation before and their first test"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " and were diagnosed with any acute exacerbation before and after their first cortisol test")))
# DF_names <- c(DF_names, ic_name)
# 
# ic_name = "No_Exacerbation"
# assign(ic_name, ICS_Tested_Compare %>% filter(Any_Acute_exacerbation == "No"))
# loginfo(str_c(nrow(get(ic_name)), " subjects had at least 1 ICS prescription, at least 1 cortisol test, but were not diagnosed with any acute exacerbation"))
# fwrite(get(ic_name), str_c(config$bonus_dir, ic_name, "_", timestamp, config$general_file_ending))
# Summary_DF <- rbind(Summary_DF,
#                     data.frame(File_Name = ic_name,
#                                nSubjects = nrow(get(ic_name)),
#                                Description = str_c("Subjects available with PPV >= 0.85,",
#                                                    " have at least 1 ICS prescription(s),",
#                                                    " had at least 1 cortisol test,",
#                                                    " but were not diagnosed with any acute exacerbation")))
# DF_names <- c(DF_names, ic_name)
# 
# 
# fwrite(Summary_DF, str_c(config$bonus_dir, "Summary_Overview", config$general_file_ending))
# 
# #(For Mengna)
# Pegasus_Ids <- fread("/pc/resta/RPDR/Pegasus/data/PEGASUS_Subject_IDs.csv") %>% pull()
# fwrite(RPDR_cleaned_data %>% filter(Biobank_Subject_ID %in% Pegasus_Ids),
#        "/pc/resta/RPDR/Pegasus/data/Pegasus_Asthma_08-31-2020_v1.csv")
# 
# ##########################
# ###### Summary Stats #####
# ##########################
# 
# Interests <- c("Age_Range", "Race", "Gender", "BMI")
# for (Category in Interests){
#   Group_Stats <- get(DF_names[1]) %>% group_by(!!as.symbol(Category)) %>%
#     summarise(!!as.symbol(str_c("Count_", DF_names[1])) := n(),
#               .groups = 'drop')
#   
#   Subcategories = DF_names[-1]
#   for (Group in Subcategories){
#     Group_Stats <- left_join(Group_Stats, get(Group) %>% group_by(!!as.symbol(Category)) %>%
#                                summarise(!!as.symbol(str_c("Count_", Group)) := n(),
#                                          .groups = 'drop'),
#                              by = Category)
#   }
#   print(Group_Stats)
# }
