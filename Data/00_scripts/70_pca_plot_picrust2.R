
#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(readxl)


#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip=1)
abundance_data  = as.data.frame(abundance_data)

#Import your metadata file, no need to filter yet
metadata <- read_xlsx("ms_metadata.xlsx")

# creating a new column with the variables (MS w/ Ashtma, MS w/o asthma, controls)
metadata <- metadata |>
  mutate(disease_var = case_when(disease=='MS' & asthma==0 ~ 'MS',
                                 disease=='MS' & asthma==1 ~ 'MS_asthma',
                                 disease=='Control' & asthma==1 ~ 'Asthma',
                                 disease=='Control' ~ 'Healthy'))



#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$disease),]

# filtering to include all 4 groups
metadata_as = metadata[!is.na(metadata$disease_var),] |>
  filter(disease_var == c('MS', 'MS_asthma', 'Healthy', 'Asthma'))

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata_as$'sample-id'
sample_names = append(sample_names, "#OTU ID")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering (based on na filtering)

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata_as = metadata_as[metadata_as$`sample-id` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                        metadata = metadata, group = "disease_var", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(`#OTU ID` %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
# 151 because we have 152 cols (inc OTU ID) in our abundance data filtered file
abundance_desc = abundance_desc[,-c(151:ncol(abundance_desc))] 

# Generate pathway PCA plot
pcaplot<- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), metadata = metadata_as, group = "disease_var")

ggsave("pcoa_plot.png", plot = pcaplot)
