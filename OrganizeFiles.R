require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(data.table) # fread, fwrite
require(lubridate) # date time conversion
require(tidyverse)
require(logging)
require(yaml)

rm(list = ls())
general_path <- "/pc/resta/RPDR/"
project_path <- str_c(general_path, "ACTH/")
setwd(project_path)
source("/pc/resta/RPDR/RPDR_General_Functions.R")

timestamp <- format(Sys.time(), "%Y-%m-%d")
config <- yaml.load_file("config.yml")
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

# Begin processing: 2 columns should be created
All_merged <- process1_biobankIDs()
# Should add 5 columns
All_merged <- process2_demographics()
# Should add 1 columns
All_merged <- process3_deidentified_ppv_only()
# Should add 9 (1 + 2*4) columns
All_merged <- process4_diagnoses(Diagnoses_Of_Interest =
                                  list("adrenal insufficiency" =
                                         c("Corticoadrenal insufficiency",
                                           "Primary adrenocortical insufficiency",
                                           "Other adrenocortical insufficiency",
                                           "Unspecified adrenocortical insufficiency")))
# Should add 57 (12 + 45) columns
All_merged <- process5_medications(Medications_Of_Interest =
                                    list("Anti-inflammatories, inhalation" =
                                           c("Beclomethasone dipropionate", "Budesonide",
                                             "Ciclesonide", "Dexamethasone",
                                             "Flunisolide", "Fluticasone",
                                             "Mometasone", "Triamcinolone"),
                                         "Antiasthma, other" = c("Fluticasone/salmeterol")),
                                   Group_Info = FALSE,
                                   merged_group_name = "Inhaled Corticosteroids")
# Should add nnnn columns
All_merged <- process6_ACTH_labs()
All_merged <- process6_ACTH_labs(strict = TRUE)
All_merged <- process7_physical()
fwrite(All_merged, output_file_name)
loginfo("Processing Complete")

##########################
###### Summary Stats #####
##########################
Summary_DF <- data.frame(File_Name = character(), nSubjects = numeric(), Description = character())

#(1) PPV 0.85
{
PPV_0.85 <- All_merged %>% filter(Asthma_0.80PPV_List_of_All_Values >= 0.85)
loginfo(str_c(nrow(PPV_0.85), " subjects available with PPV >= 0.85"))
fwrite(PPV_0.85, str_c(config$data_dir, "PPV_0.85_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = PPV_0.85, additional_header = "PPV_0.85_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "PPV_0.85",
                               nSubjects = nrow(PPV_0.85),
                               Description = "Subjects available with PPV >= 0.85"))
}

#(2) ICS use
{
ICS <- PPV_0.85 %>% filter(grepl("Yes", Any_Inhaled_Corticosteroids))
loginfo(str_c(nrow(ICS), " subjects have at least 1 ICS prescription"))
fwrite(ICS, str_c(config$data_dir, "ICS_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS, additional_header = "ICS_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS",
                               nSubjects = nrow(ICS),
                               Description = str_c("Subjects available with PPV >= 0.85",
                                                   " and have at least 1 ICS prescription")))

ICS_4 <- ICS %>% filter(Any_Inhaled_Corticosteroids_total_dates >= 4)
loginfo(str_c(nrow(ICS_4), " subjects have at least 4 ICS prescriptions"))
fwrite(ICS_4, str_c(config$data_dir, "ICS_4_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4, additional_header = "ICS_4_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_4",
                               nSubjects = nrow(ICS_4),
                               Description = str_c("Subjects available with PPV >= 0.85",
                                                   " and have at least 4 ICS prescriptions")))

ICS_10 <- ICS %>% filter(Any_Inhaled_Corticosteroids_total_dates >= 10)
loginfo(str_c(nrow(ICS_10), " subjects have at least 10 ICS prescriptions"))
fwrite(ICS_10, str_c(config$data_dir, "ICS_10_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10, additional_header = "ICS_10_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_10",
                               nSubjects = nrow(ICS_10),
                               Description = str_c("Subjects available with PPV >= 0.85",
                                                   " and have at least 10 ICS prescriptions")))
}

#(3a) ICS use and cortisol test
{
ICS_Tested <- ICS %>% filter(All_Cortisol_nTotalDates >= 1)
loginfo(str_c(nrow(ICS_Tested), " subjects have at least 1 ICS prescription and had at least 1 cortisol test"))
fwrite(ICS_Tested, str_c(config$data_dir, "ICS_Tested_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested, additional_header = "ICS_Tested_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_Tested",
                               nSubjects = nrow(ICS_Tested),
                               Description = str_c("Subjects available with PPV >= 0.85,",
                                                   " have at least 1 ICS prescription,",
                                                   " and have at least one cortisol test")))

ICS_4_Tested <- ICS_4 %>% filter(All_Cortisol_nTotalDates >= 1)
loginfo(str_c(nrow(ICS_4_Tested), " subjects have at least 4 ICS prescriptions and had at least 1 cortisol test"))
fwrite(ICS_4_Tested, str_c(config$data_dir, "ICS_4_Tested_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested, additional_header = "ICS_4_Tested_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_4_Tested",
                               nSubjects = nrow(ICS_4_Tested),
                               Description = str_c("Subjects available with PPV >= 0.85,",
                                                   " have at least 4 ICS prescriptions,",
                                                   " and have at least one cortisol test")))

ICS_10_Tested <- ICS_10 %>% filter(All_Cortisol_nTotalDates >= 1)
loginfo(str_c(nrow(ICS_10_Tested), " subjects have at least 10 ICS prescriptions and had at least 1 cortisol test"))
fwrite(ICS_10_Tested, str_c(config$data_dir, "ICS_10_Tested_", timestamp, config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested, additional_header = "ICS_10_Tested_", return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_10_Tested",
                               nSubjects = nrow(ICS_10_Tested),
                               Description = str_c("Subjects available with PPV >= 0.85,",
                                                   " have at least 10 ICS prescriptions,",
                                                   " and have at least one cortisol test")))
}

#(3b) ICS use and cortisol test with cutoffs
{
  ##(3b.1) Cutoff 2
  {
    Cutoff = 2
    ICS_Tested_Cutoff_2 <- ICS_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_Tested_Cutoff_2), " subjects have at least 1 ICS prescription,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_Tested_Cutoff_2, str_c(config$data_dir, "ICS_Tested_Cutoff_", Cutoff, "_", timestamp,
                                      config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested_Cutoff_2,
                               additional_header = str_c("ICS_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_Tested_Cutoff_2",
                                   nSubjects = nrow(ICS_Tested_Cutoff_2),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 1 ICS prescription,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    
    ICS_4_Tested_Cutoff_2 <- ICS_4_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_4_Tested_Cutoff_2), " subjects have at least 4 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_4_Tested_Cutoff_2, str_c(config$data_dir, "ICS_4_Tested_Cutoff_", Cutoff, "_", timestamp,
                                        config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested_Cutoff_2,
                               additional_header = str_c("ICS_4_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_4_Tested_Cutoff_2",
                                   nSubjects = nrow(ICS_4_Tested_Cutoff_2),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 4 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_10_Tested_Cutoff_2 <- ICS_10_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_10_Tested_Cutoff_2), " subjects have at least 10 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_10_Tested_Cutoff_2, str_c(config$data_dir, "ICS_10_Tested_Cutoff_", Cutoff, "_", timestamp,
                                         config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested_Cutoff_2,
                               additional_header = str_c("ICS_10_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_10_Tested_Cutoff_2",
                                   nSubjects = nrow(ICS_10_Tested_Cutoff_2),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 10 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
  }
  ##(3b.2) Cutoff 5
  {
    Cutoff = 5
    ICS_Tested_Cutoff_5 <- ICS_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_Tested_Cutoff_5), " subjects have at least 1 ICS prescription,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_Tested_Cutoff_5, str_c(config$data_dir, "ICS_Tested_Cutoff_", Cutoff, "_", timestamp,
                                      config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested_Cutoff_5,
                               additional_header = str_c("ICS_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_Tested_Cutoff_5",
                                   nSubjects = nrow(ICS_Tested_Cutoff_5),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 1 ICS prescription,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_4_Tested_Cutoff_5 <- ICS_4_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_4_Tested_Cutoff_5), " subjects have at least 4 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_4_Tested_Cutoff_5, str_c(config$data_dir, "ICS_4_Tested_Cutoff_", Cutoff, "_", timestamp,
                                        config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested_Cutoff_5,
                               additional_header = str_c("ICS_4_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_4_Tested_Cutoff_5",
                                   nSubjects = nrow(ICS_4_Tested_Cutoff_5),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 4 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_10_Tested_Cutoff_5 <- ICS_10_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_10_Tested_Cutoff_5), " subjects have at least 10 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_10_Tested_Cutoff_5, str_c(config$data_dir, "ICS_10_Tested_Cutoff_", Cutoff, "_", timestamp,
                                         config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested_Cutoff_5,
                               additional_header = str_c("ICS_10_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_10_Tested_Cutoff_5",
                                   nSubjects = nrow(ICS_10_Tested_Cutoff_5),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 10 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
  }
  ##(3b.3) Cutoff 10
  {
    Cutoff = 10
    ICS_Tested_Cutoff_10 <- ICS_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_Tested_Cutoff_10), " subjects have at least 1 ICS prescription,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_Tested_Cutoff_10, str_c(config$data_dir, "ICS_Tested_Cutoff_", Cutoff, "_", timestamp,
                                       config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested_Cutoff_10,
                               additional_header = str_c("ICS_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_Tested_Cutoff_10",
                                   nSubjects = nrow(ICS_Tested_Cutoff_10),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 1 ICS prescription,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_4_Tested_Cutoff_10 <- ICS_4_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_4_Tested_Cutoff_10), " subjects have at least 4 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_4_Tested_Cutoff_10, str_c(config$data_dir, "ICS_4_Tested_Cutoff_", Cutoff, "_", timestamp,
                                         config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested_Cutoff_10,
                               additional_header = str_c("ICS_4_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_4_Tested_Cutoff_10",
                                   nSubjects = nrow(ICS_4_Tested_Cutoff_10),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 4 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_10_Tested_Cutoff_10 <- ICS_10_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_10_Tested_Cutoff_10), " subjects have at least 10 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_10_Tested_Cutoff_10, str_c(config$data_dir, "ICS_10_Tested_Cutoff_", Cutoff, "_", timestamp,
                                          config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested_Cutoff_10,
                               additional_header = str_c("ICS_10_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_10_Tested_Cutoff_10",
                                   nSubjects = nrow(ICS_10_Tested_Cutoff_10),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 10 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
  }
  ##(3b.4) Cutoff 15
  {
    Cutoff = 15
    ICS_Tested_Cutoff_15 <- ICS_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_Tested_Cutoff_15), " subjects have at least 1 ICS prescription,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_Tested_Cutoff_15, str_c(config$data_dir, "ICS_Tested_Cutoff_", Cutoff, "_", timestamp,
                                       config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested_Cutoff_15,
                               additional_header = str_c("ICS_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_Tested_Cutoff_15",
                                   nSubjects = nrow(ICS_Tested_Cutoff_15),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 1 ICS prescription,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_4_Tested_Cutoff_15 <- ICS_4_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_4_Tested_Cutoff_15), " subjects have at least 4 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_4_Tested_Cutoff_15, str_c(config$data_dir, "ICS_4_Tested_Cutoff_", Cutoff, "_", timestamp,
                                         config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested_Cutoff_15,
                               additional_header = str_c("ICS_4_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_4_Tested_Cutoff_15",
                                   nSubjects = nrow(ICS_4_Tested_Cutoff_15),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 4 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_10_Tested_Cutoff_15 <- ICS_10_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_10_Tested_Cutoff_15), " subjects have at least 10 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_10_Tested_Cutoff_15, str_c(config$data_dir, "ICS_10_Tested_Cutoff_", Cutoff, "_", timestamp,
                                          config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested_Cutoff_15,
                               additional_header = str_c("ICS_10_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_10_Tested_Cutoff_15",
                                   nSubjects = nrow(ICS_10_Tested_Cutoff_15),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 10 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
  }
  ##(3b.5) Cutoff 200
  {
    Cutoff = 20
    ICS_Tested_Cutoff_20 <- ICS_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_Tested_Cutoff_20), " subjects have at least 1 ICS prescription,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_Tested_Cutoff_20, str_c(config$data_dir, "ICS_Tested_Cutoff_", Cutoff, "_", timestamp,
                                       config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested_Cutoff_20,
                               additional_header = str_c("ICS_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_Tested_Cutoff_20",
                                   nSubjects = nrow(ICS_Tested_Cutoff_20),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 1 ICS prescription,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_4_Tested_Cutoff_20 <- ICS_4_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_4_Tested_Cutoff_20), " subjects have at least 4 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_4_Tested_Cutoff_20, str_c(config$data_dir, "ICS_4_Tested_Cutoff_", Cutoff, "_", timestamp,
                                         config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested_Cutoff_20,
                               additional_header = str_c("ICS_4_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_4_Tested_Cutoff_20",
                                   nSubjects = nrow(ICS_4_Tested_Cutoff_20),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 4 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
    
    ICS_10_Tested_Cutoff_20 <- ICS_10_Tested %>% filter(All_Cortisol_Overall_Median_Result <= Cutoff)
    loginfo(str_c(nrow(ICS_10_Tested_Cutoff_20), " subjects have at least 10 ICS prescriptions,",
                  "had at least 1 cortisol test, and the median value was no greater than ", Cutoff))
    fwrite(ICS_10_Tested_Cutoff_20, str_c(config$data_dir, "ICS_10_Tested_Cutoff_", Cutoff, "_", timestamp,
                                          config$general_file_ending))
    create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested_Cutoff_20,
                               additional_header = str_c("ICS_10_Tested_Cutoff_", Cutoff, "_"),
                               return.cholesterol = FALSE)
    Summary_DF <- rbind(Summary_DF,
                        data.frame(File_Name = "ICS_10_Tested_Cutoff_20",
                                   nSubjects = nrow(ICS_10_Tested_Cutoff_20),
                                   Description = str_c("Subjects available with PPV >= 0.85,",
                                                       " have at least 10 ICS prescriptions,",
                                                       " and have at least one cortisol test",
                                                       " where median is not greater than ", Cutoff)))
  }
}

#(4) ICS use and cortisol test and adrenal diagnosis
{
ICS_Tested_Diagnosis <- ICS_Tested %>% filter(Any_adrenal_insufficiency_diagnosis == "Yes")
loginfo(str_c(nrow(ICS_Tested_Diagnosis), " subjects have at least 1 ICS prescription,",
              " had at least 1 cortisol test, and at least one adrenal insufficiency diagnosis"))
fwrite(ICS_Tested_Diagnosis, str_c(config$data_dir, "ICS_Tested_Diagnosis_", timestamp,
                                   config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_Tested_Diagnosis,
                           additional_header = "ICS_Tested_Diagnosis_",
                           return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_Tested_Diagnosis",
                               nSubjects = nrow(ICS_Tested_Diagnosis),
                               Description = str_c("Subjects available with PPV >= 0.85,",
                                                   " have at least 1 ICS prescription,",
                                                   " have at least one cortisol test,",
                                                   " and have at least 1 adrenal insufficency diagnosis")))

ICS_4_Tested_Diagnosis <- ICS_4_Tested %>% filter(Any_adrenal_insufficiency_diagnosis == "Yes")
loginfo(str_c(nrow(ICS_4_Tested_Diagnosis), " subjects have at least 4 ICS prescriptions,",
              " had at least 1 cortisol test, and at least one adrenal insufficiency diagnosis"))
fwrite(ICS_4_Tested_Diagnosis, str_c(config$data_dir, "ICS_4_Tested_Diagnosis_", timestamp,
                                     config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_4_Tested_Diagnosis,
                           additional_header = "ICS_4_Tested_Diagnosis_",
                           return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_4_Tested_Diagnosis",
                               nSubjects = nrow(ICS_4_Tested_Diagnosis),
                               Description = str_c("Subjects available with PPV >= 0.85,",
                                                   " have at least 4 ICS prescriptions,",
                                                   " have at least one cortisol test,",
                                                   " and have at least 1 adrenal insufficency diagnosis")))

ICS_10_Tested_Diagnosis <- ICS_10_Tested %>% filter(Any_adrenal_insufficiency_diagnosis == "Yes")
loginfo(str_c(nrow(ICS_10_Tested_Diagnosis), " subjects have at least 10 ICS prescriptions,",
              " had at least 1 cortisol test, and at least one adrenal insufficiency diagnosis"))
fwrite(ICS_10_Tested_Diagnosis, str_c(config$data_dir, "ICS_10_Tested_Diagnosis_", timestamp,
                                      config$general_file_ending))
create_summary_stats_files(use.TS = FALSE, Cleaned_DF = ICS_10_Tested_Diagnosis,
                           additional_header = "ICS_10_Tested_Diagnosis_",
                           return.cholesterol = FALSE)
Summary_DF <- rbind(Summary_DF,
                    data.frame(File_Name = "ICS_10_Tested_Diagnosis",
                               nSubjects = nrow(ICS_10_Tested_Diagnosis),
                               Description = str_c("Subjects available with PPV >= 0.85,",
                                                   " have at least 10 ICS prescriptions,",
                                                   " have at least one cortisol test,",
                                                   " and have at least 1 adrenal insufficency diagnosis")))
}

#(5) ICS changes before and after cortisol test
names(ICS)


fwrite(Summary_DF, str_c(config$data_dir, "Summary_Overview", config$general_file_ending))
