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
config <- yaml.load_file("configACTHAdrenalQuery.yml")
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
source(str_c(config$rpdr_src_dir, "Clean_Health_History_File.R")) # process_physical
source(str_c(config$rpdr_src_dir, "Clean_Diagnoses_File.R")) # process_diagnoses
source(str_c(config$rpdr_src_dir, "Clean_Medications_File.R")) # process_medications
source(str_c(config$rpdr_src_dir, "Clean_Labs_File.R")) # process_labs

# Begin processing and cleaning of RPDR files: 10 columns should be created
All_merged <- start_processing()

# Should add 8 columns -> 18 columns
All_merged <- process_deidentified(clean.list = TRUE)

# Should add 4 columns -> 22 columns
All_merged <- process_plasma(date_cutoff = "2020-06-30", reduce_to_plasma_only = FALSE)

# Get file of Genotype Ids 2 more columns -> 24
Genotype_File <- fread(config$genotype_file_name)
Genotype_File <- Genotype_File %>% rename("Biobank_Subject_ID" = SubjectID)
All_merged <- left_join(All_merged, Genotype_File) %>%
  mutate(Has_Genotype_Info_MegE = ifelse(is.na(SentrixID), "0", "1")) %>%
  select(EMPI:All_Collection_Dates, Has_Genotype_Info_MegE, SentrixID, everything())
rm(Genotype_File)

# Should add 42 columns -> 66 columns
All_merged <- process_physical(Return_Influenza = FALSE,
                               Return_Blood_Pressure = TRUE)

# Should add 26 (6 + 4*5) columns -> 92 columns
All_merged <- process_diagnoses(Diagnoses_Of_Interest =
                                  list("Adrenal insufficiency" =
                                         c("Corticoadrenal insufficiency",
                                           "Primary adrenocortical insufficiency",
                                           "Other adrenocortical insufficiency",
                                           "Unspecified adrenocortical insufficiency")))

# Should add 36 (5*6 + 6) columns -> 128 columns
All_merged <- process_diagnoses(Diagnoses_Of_Interest =
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
                                Merge_Group_Info_Name = "Asthma")

# Should add 12 (2*6) columns -> 140 columns
All_merged <- process_diagnoses(Diagnoses_Of_Interest =
                                  list("Bronchiestasis" = c("^Bronchiectasis"),
                                       "Other chronic obstructive pulmonary disease" =
                                         c("Chronic airway",
                                           "Chronic bronchitis-",
                                           "Chronic obstructive asthma",
                                           "Chronic obstructive lung",
                                           "Chronic obstructive pulmonary disease.[^\\(]",
                                           "Obstructive chronic bronchitis")),
                                Exact = FALSE,
                                Individual_Info = FALSE)

# Should add 48 (7*6 + 6) columns -> 188 columns
All_merged <- process_diagnoses(Diagnoses_Of_Interest =
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
                                Merge_Group_Info_Name = "Acute exacerbation")

# Should add 100 (9 + 13*7) columns -> 288 columns
All_merged <- process_medications(Medications_Of_Interest =
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
                                  merged_group_name = "Inhaled Corticosteroids")

# Should add 142 (9 + 19*7) columns -> 430 columns
All_merged <- process_medications(Medications_Of_Interest =
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
                                  merged_group_name = "Not Inhaled Corticosteroids")

# Should add 210 columns -> 640 columns
All_merged <- process_ACTH_labs()

# Relevant previous information ends here

# New changes for comparison
# Should add 56 -> 696 columns
All_merged <- process_diagnoses(Diagnoses_Of_Interest =
                                  list("Fatigue" =
                                         c("[Mm]alaise and fatigue",
                                           "Chronic fatigue",
                                           "Fatigue-LMR",
                                           "Fatigue-Oncall",
                                           "Neoplastic.*fatigue",
                                           "Other.*fatigue"),
                                       "Anemia" =
                                         c("[Aa]cquired hemolytic anemia",
                                           "Acute posthemorrhagic anemia",
                                           "Anemia associated with other specified nutritional deficiency",
                                           "Anemia(s)* due to",
                                           "Anemia in",
                                           "Anemia of (other ){0,1}chronic",
                                           "Anemia-",
                                           "Anemia, unspecified",
                                           "Antineoplastic chemotherapy induced anemia",
                                           "[Aa]plastic anemia",
                                           "[Aa]utoimmune hemolytic anemia",
                                           "[Bb]12 deficiency anemia",
                                           "Congenital dyserythropoietic anemia",
                                           "Drug-induced.*anemia",
                                           "[Ff]olate.deficiency anemia",
                                           "Hemolytic anemia",
                                           "Hereditary.*anemia",
                                           "Iron deficiency anemia",
                                           "[Mm]icrocytic anemia",
                                           "[Nn]utritional anemia",
                                           "Other .* anemia",
                                           "Pernicious anemia",
                                           "Protein deficiency anemia",
                                           "[Ss]icle.cell anemia",
                                           "Sideroblastic anemia",
                                           "Unspecified deficiency anemia"),
                                       "Weight Loss" = 
                                         c("[Ll]oss of weight",
                                           "[Ww]eight loss",
                                           "Poor weight gain"),
                                       "Hyperpigmentation" = 
                                         c("[Dd]isorder(s)* of pigmentation",
                                           "[Dd]yschromia",
                                           "Hyperpigmentation-Oncall",
                                           "hyperpigmentation",
                                           "Other disorders of skin and subcutaneous tissue",
                                           "Other melanin"),
                                       "Lightheadedness" = 
                                         c("Lightheadedness")),
                                Individual_Info = FALSE,
                                Exact = FALSE,
                                Merge_Group_Info_Name = "Adrenal Fatigue")

# Filter data so that it's only subjects with any ACTH test
Cortisol <- All_merged %>% filter(!is.na(All_Cortisol_nTotalDates))
# Add 2 columns -> 688
Cortisol <- Cortisol %>%
  mutate(MostRecentMinus1Year = All_Cortisol_Date_Last - years(1),
         MostRecentMinus5Years = All_Cortisol_Date_Last - years(5))
# 
# # Should add 28 columns -> 49 columns
# All_merged <- process_physical_set_range(Return_Influenza = FALSE,
#                                Return_Blood_Pressure = TRUE,
#                                min_dates = "MostRecentMinus1Year",
#                                max_dates = "All_Cortisol_Date_Last",
#                                Range_Name = "OneYearBeforeTest")
# 
# # Should add 28 columns -> 49 columns
# All_merged <- process_physical_set_range(Return_Influenza = FALSE,
#                                Return_Blood_Pressure = TRUE,
#                                min_dates = "MostRecentMinus5Years",
#                                max_dates = "All_Cortisol_Date_Last",
#                                Range_Name = "FiveYearsBeforeTest")

# Should add 26 (6 + 4*5) columns -> 712 columns
Cortisol <- process_diagnoses_set_range(DF_to_fill = Cortisol,
                                        Diagnoses_Of_Interest =
                                          list("Adrenal insufficiency" =
                                                 c("Corticoadrenal insufficiency",
                                                   "Primary adrenocortical insufficiency",
                                                   "Other adrenocortical insufficiency",
                                                   "Unspecified adrenocortical insufficiency")),
                                        min_dates = "MostRecentMinus1Year",
                                        max_dates = "All_Cortisol_Date_Last",
                                        Range_Name = "OneYearBeforeTest")

# Should add 26 (6 + 4*5) columns -> 712 columns
Cortisol <- process_diagnoses_set_range(DF_to_fill = Cortisol,
                                        Diagnoses_Of_Interest =
                                          list("Adrenal insufficiency" =
                                                 c("Corticoadrenal insufficiency",
                                                   "Primary adrenocortical insufficiency",
                                                   "Other adrenocortical insufficiency",
                                                   "Unspecified adrenocortical insufficiency")),
                                        min_dates = "MostRecentMinus5Years",
                                        max_dates = "All_Cortisol_Date_Last",
                                        Range_Name = "FiveYearsBeforeTest")

# New changes for comparison
# Should add 26 (6 + 4*5) columns -> 75 columns
Cortisol <- process_diagnoses_set_range(DF_to_fill = Cortisol,
                                          Diagnoses_Of_Interest =
                                            list("Fatigue" =
                                                   c("[Mm]alaise and fatigue",
                                                     "Chronic fatigue",
                                                     "Fatigue-LMR",
                                                     "Fatigue-Oncall",
                                                     "Neoplastic.*fatigue",
                                                     "Other.*fatigue"),
                                                 "Anemia" =
                                                   c("[Aa]cquired hemolytic anemia",
                                                     "Acute posthemorrhagic anemia",
                                                     "Anemia associated with other specified nutritional deficiency",
                                                     "Anemia(s)* due to",
                                                     "Anemia in",
                                                     "Anemia of (other ){0,1}chronic",
                                                     "Anemia-",
                                                     "Anemia, unspecified",
                                                     "Antineoplastic chemotherapy induced anemia",
                                                     "[Aa]plastic anemia",
                                                     "[Aa]utoimmune hemolytic anemia",
                                                     "[Bb]12 deficiency anemia",
                                                     "Congenital dyserythropoietic anemia",
                                                     "Drug-induced.*anemia",
                                                     "[Ff]olate.deficiency anemia",
                                                     "Hemolytic anemia",
                                                     "Hereditary.*anemia",
                                                     "Iron deficiency anemia",
                                                     "[Mm]icrocytic anemia",
                                                     "[Nn]utritional anemia",
                                                     "Other .* anemia",
                                                     "Pernicious anemia",
                                                     "Protein deficiency anemia",
                                                     "[Ss]icle.cell anemia",
                                                     "Sideroblastic anemia",
                                                     "Unspecified deficiency anemia"),
                                                 "Weight Loss" = 
                                                   c("[Ll]oss of weight",
                                                     "[Ww]eight loss",
                                                     "Poor weight gain"),
                                                 "Hyperpigmentation" = 
                                                   c("[Dd]isorder(s)* of pigmentation",
                                                     "[Dd]yschromia",
                                                     "Hyperpigmentation-Oncall",
                                                     "hyperpigmentation",
                                                     "Other disorders of skin and subcutaneous tissue",
                                                     "Other melanin"),
                                                 "Lightheadedness" = 
                                                   c("Lightheadedness")),
                                          Individual_Info = FALSE,
                                          Exact = FALSE,
                                          Merge_Group_Info_Name = "Adrenal Fatigue",
                                          min_dates = "MostRecentMinus1Year",
                                          max_dates = "All_Cortisol_Date_Last",
                                          Range_Name = "OneYearBeforeTest")

# New changes for comparison
# Should add 26 (6 + 4*5) columns -> 75 columns
Cortisol <- process_diagnoses_set_range(DF_to_fill = Cortisol,
                                          Diagnoses_Of_Interest =
                                            list("Fatigue" =
                                                   c("[Mm]alaise and fatigue",
                                                     "Chronic fatigue",
                                                     "Fatigue-LMR",
                                                     "Fatigue-Oncall",
                                                     "Neoplastic.*fatigue",
                                                     "Other.*fatigue"),
                                                 "Anemia" =
                                                   c("[Aa]cquired hemolytic anemia",
                                                     "Acute posthemorrhagic anemia",
                                                     "Anemia associated with other specified nutritional deficiency",
                                                     "Anemia(s)* due to",
                                                     "Anemia in",
                                                     "Anemia of (other ){0,1}chronic",
                                                     "Anemia-",
                                                     "Anemia, unspecified",
                                                     "Antineoplastic chemotherapy induced anemia",
                                                     "[Aa]plastic anemia",
                                                     "[Aa]utoimmune hemolytic anemia",
                                                     "[Bb]12 deficiency anemia",
                                                     "Congenital dyserythropoietic anemia",
                                                     "Drug-induced.*anemia",
                                                     "[Ff]olate.deficiency anemia",
                                                     "Hemolytic anemia",
                                                     "Hereditary.*anemia",
                                                     "Iron deficiency anemia",
                                                     "[Mm]icrocytic anemia",
                                                     "[Nn]utritional anemia",
                                                     "Other .* anemia",
                                                     "Pernicious anemia",
                                                     "Protein deficiency anemia",
                                                     "[Ss]icle.cell anemia",
                                                     "Sideroblastic anemia",
                                                     "Unspecified deficiency anemia"),
                                                 "Weight Loss" = 
                                                   c("[Ll]oss of weight",
                                                     "[Ww]eight loss",
                                                     "Poor weight gain"),
                                                 "Hyperpigmentation" = 
                                                   c("[Dd]isorder(s)* of pigmentation",
                                                     "[Dd]yschromia",
                                                     "Hyperpigmentation-Oncall",
                                                     "hyperpigmentation",
                                                     "Other disorders of skin and subcutaneous tissue",
                                                     "Other melanin"),
                                                 "Lightheadedness" = 
                                                   c("Lightheadedness")),
                                          Individual_Info = FALSE,
                                          Exact = FALSE,
                                          Merge_Group_Info_Name = "Adrenal Fatigue",
                                          min_dates = "MostRecentMinus5Years",
                                          max_dates = "All_Cortisol_Date_Last",
                                          Range_Name = "FiveYearsBeforeTest")

All_merged <- left_join(All_merged, Cortisol)

fwrite(All_merged, output_file_name)
loginfo("Processing Complete")
