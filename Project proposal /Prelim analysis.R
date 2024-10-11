library(tidyverse)
library(readxl)
#Loading data
data <- "ms_metadata.xlsx"
ms_data <- read_excel(data)
ms_data

#Renaming sample ID column
colnames(ms_data)[1] <- "sample_id" 

#Selecting for MS and asthma
selected <- select(ms_data, sample_id, disease_course, asthma)
selected

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
