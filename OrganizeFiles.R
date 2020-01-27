require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(data.table) # fread, fwrite
require(lubridate) # date time conversion
require(tidyverse)
require(logging)
require(yaml)

setwd("/pc/resta/RPDR/ACTH")
# getwd()
# list.files()
rm(list = ls())
source("/pc/resta/RPDR/RPDR_General_Functions.R")

timestamp <- format(Sys.time(), "%Y-%m-%d")
config <- yaml.load_file("config.yml")
output_file_name = str_c(config$output_file_header, timestamp, config$general_file_ending)
if (config$create_intermediates){
  if (!dir.exists(config$select_files_dir)){
    dir.create(config$select_files_dir, recursive = TRUE)
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
# Should add 2 columns
All_merged <- process3_deidentified()
# Should add 9 (1 + 2*4) columns
All_merged <- process4_diagnoses(Diagnoses_Of_Interest =
                                  list("adrenal insufficiency" =
                                         c("Corticoadrenal insufficiency",
                                           "Primary adrenocortical insufficiency",
                                           "Other adrenocortical insufficiency",
                                           "Unspecified adrenocortical insufficiency")))
# Should add 57 (12 + 45) columns
All_merged <- process5_medications(Medications_Of_Interest =
                                    list("Inhaled Corticosteroids" =
                                           c("Beclomethasone dipropionate", "Budesonide",
                                             "Ciclesonide", "Dexamethasone",
                                             "Flunisolide", "Fluticasone",
                                             "Fluticasone/salmeterol", "Mometasone",
                                             "Triamcinolone")))
# Should add nnnn columns
All_merged <- process6_ACTH_labs()
All_merged <- process6_ACTH_labs(strict = TRUE)
fwrite(All_merged, output_file_name)
loginfo("Processing Complete")
