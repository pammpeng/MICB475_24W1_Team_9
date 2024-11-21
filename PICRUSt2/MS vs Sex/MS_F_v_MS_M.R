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

# MS_male vs MS_female
#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip=1)
abundance_data  =as.data.frame(abundance_data) %>% 
  rename('#OTU ID' = 'pathway')

#Import your metadata file, no need to filter yet
metadata <- read_delim("ms_metadata_final_2.tsv")

# creating a new column with the variables 
metadata <- metadata |>
  mutate(disease = case_when(disease_course=='RRMS' & sex=='F' ~ 'MS_F',
                             disease_course=='RRMS' & sex=='M' ~ 'MS_M',
                             disease_course=='Control' & sex=='F' ~ 'Control_F',
                             disease_course=='Control' & sex=='M' ~ 'Control_M'))

#Example Looking at subject number
#If you have multiple variants, filter your metadata to include only 2 at a time

metadata <- metadata %>% 
  filter(disease == c('MS_M', 'MS_F'))

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$disease),]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata$'sample-id'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata = metadata[metadata$`sample-id` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                        metadata = metadata, group = "disease", daa_method = "DESeq2")

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
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 38 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(37:ncol(abundance_desc))] 

# Generate a heatmap
heatmap <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata, group = "disease")
heatmap2<- heatmap + theme(axis.text.y = element_text(size = 6)) 
ggsave("heatmap_msvssex.png", plot = heatmap2, width = 10, height = 6)
# Generate pathway PCA plot
pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadata, group = "disease")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")


# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "disease")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05 & log2FoldChange > -1.5 & log2FoldChange < 1.5)
# You can also filter by Log2fold change

sig_res <- sig_res[order(sig_res$log2FoldChange),]
sig_res2<-ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")
ggsave("sigres_MSvssex.png", plot = sig_res2, width = 10, height = 6)
