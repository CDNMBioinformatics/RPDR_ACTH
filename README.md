# RPDR_Priya

This repository contains the code for parsing the Asthma and ACTH data from Biobank and RPDR. This may include future analysis code but it depends on where the project goes.

To find identified data, go to /pc/resta/RPDR/Project_For_Priya/data . Note only authorized users have access to this directory.

Relevant data for this project which is not checked into git:
- BiobankPortal_mns43_2019-09-18-150014.csv # Includes Asthma and Nonasthma variables
- mns43_100719120614239897_Med.txt # Includes medication information
- mns43_100719120614239897_Dem.txt # Includes demographic information
- mns43_100719120614239897_Dia.txt # Includes diagnostic information
- mns43_100719120614239897_Lab.txt # Includes requested lab information

Checked in files and directories:
- OrganizeFiles.R # Read in each of the RPDR files and Biobank files to clean and parse for relevant study information
- Reports/ # pdf and text files listed what was queried and delivered

NOTE: This project is (or will be) a submodule of resta/RPDR which includes the medication mapping function
