library(tidyverse)
library(readxl)

# renaming sample-id
colnames(ms_metadata)[1] <- 'sample_id'

cols_kept <- c('sample_id', 'site_x', 'sex', 'age',
               'bmi', 'disease_course', 'treatment_status', 'treatments',
               'administration', 'MSSS', 'dmt','birth_method', 'breastfeeding',
               'allergies','asthma', 'ms_family', 'nsaids', 'smoke', 'education',
               'occupation', 'vitamin D (IU)')

filtered_metadata <- ms_metadata|>
  select(all_of(cols_kept)) |>
  distinct(sample_id, .keep_all = TRUE) # removing duplicates

# tb for MS

ms_tb <- filtered_metadata |>
  mutate(disease_course = ifelse(disease_course!='Control', 'MS', disease_course))

# Fisher's or Chi-squared test test -------------------------------------

# sex
fisher_sex <- fisher.test(table(ms_tb$sex, ms_tb$disease_course))

# education
chisq_ed <- chisq.test(table(ms_tb$education, ms_tb$disease_course))

# birth method
chisq_bm <- chisq.test(table(ms_tb$birth_method, ms_tb$disease_course))


# filtering for non-control 
nc_tb <- filtered_metadata|>
  filter(disease_course != 'Control')
