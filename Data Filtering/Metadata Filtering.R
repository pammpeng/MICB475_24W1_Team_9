library(tidyverse)
library(readxl)
library(dplyr)
#Loading data
data <- "ms_metadata.xlsx"
ms_data <- read_excel(data)
ms_data

#### Preliminary Filtering to Determine Sample Size ####

#filtering for healthy and asthma

#Selecting for MS and asthma
selected <- select(ms_data, "sample-id", "disease_course", "asthma")
selected

#filtering for healthy and asthma
filtered_asthma <- filter(selected, disease_course == "Control", asthma == 1)
filtered_asthma

#Filtering for PPMS 
filtered_PPMS <- filter(selected, disease_course == "PPMS")
filtered_PPMS

#Filtering for asthma within PPMS 
filtered_PPMS_asthma <- filter(filtered_PPMS, asthma == 1)
filtered_PPMS_asthma

#Filtering for RRMS
rrms_filtered <- filter(selected, disease_course == "RRMS")
rrms_filtered

#Filtering for asthma within RRMS
rrms_asthma_filtered <- filter(rrms_filtered, asthma == 1)
rrms_asthma_filtered

#Removing duplicated sample IDs
rrms_asthma_filtered_unique <- unique(rrms_asthma_filtered)
rrms_asthma_filtered_unique

#Counting the number of rows from both RRMS and PPMS + asthma
nrow(rrms_asthma_filtered_unique)
nrow(filtered_PPMS_asthma)
nrow(filtered_asthma)

#### Removing Unwanted Variables and Selecting Samples for Analysis ####
# remove duplicates from data 
ms_no_dbl <- unique(ms_data)

#remove unwanted variables 
ms_clean <- select(ms_no_dbl, "sample-id","sex","age",
                   "disease_course","treatment_status",
                   "asthma","smoke")

#create RRMS and asthma group 
rrms_asthma <- filter(ms_clean, disease_course =="RRMS", asthma == 1)

#create RRMS, no asthma group 
set.seed(9)
rrms_no_asthma <- filter(ms_clean, disease_course =="RRMS", asthma ==0) %>%
  sample_n(38)

#create healthy, asthma group 
healthy_asthma <- filter(ms_clean, disease_course =="Control", asthma ==1) %>%
  sample_n(38)

#create healthy, no asthma group 
healthy_no_asthma <- filter(ms_clean, disease_course =="Control", asthma ==0) %>%
  sample_n(38)

#combine all 4 groups into 1 data frame
ms_final_dat <- full_join(rrms_asthma, rrms_no_asthma) %>%
  full_join(healthy_asthma) %>%
  full_join(healthy_no_asthma)
  
#export the final data frame as a tsv file and a csv file
write_delim(ms_final_dat, "ms_metadata_filtered.tsv", delim = "\t")

