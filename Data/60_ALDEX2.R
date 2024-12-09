# Load libraries and data  

install.packages("randomcoloR")
install.packages("ggpubr")
library(tidyverse) # For all your basic data wrangling and plotting needs.
library(phyloseq) # Indispensable package for microbiome analyses. 
library(ggpubr)
library(randomcoloR)
library(purrr)

# Create a folder for plot .jpegs, if it doesn't already exist
if (!file.exists('Plots')){ dir.create('Plots')}

#Setting the seed to ensure reproducible analyses
set.seed(711) 

# Initial Data Processing
#Loading the phyloseq object
load("MS_phyloseq.Rdata")

# Extract metadata (sample data) from phyloseq object
metadata <- sample_data(MS_phyloseq)

# Create the new column based on conditions in the 'disease_course' and 'asthma' columns
metadata$MS_asthma <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 1, "Yes", "No")

metadata$MS_asthma_vs_asthma <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 1, "MS_asthma", 
                                       ifelse(metadata$asthma == 1, "asthma", "No"))

metadata$MS_asthma_vs_MS <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 1, "MS_asthma", 
                                   ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 0, "MS", NA))

metadata$MS_vs_asthma <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 0, "MS", 
                                ifelse(metadata$disease_course == "Control" & metadata$asthma == 1, "asthma", NA))

# Check the first few rows of the updated metadata
head(metadata)

# Add the new column back to the sample data of the phyloseq object
sample_data(MS_phyloseq) <- metadata

# Check if the new column is added to the phyloseq object
sample_data(MS_phyloseq)

# Let's see how our sequencing depth looks
hist(sample_sums(MS_phyloseq),breaks = 30) 
table(below_1000 = sample_sums(MS_phyloseq)<=1000) 
#  expt = MS_phyloseq@sam_data$disease_course) # shows us how many samples are below 1000
#All 151 samples have more than 1000 reads --> no need to further prune
#75 samples in the Control group (no MS) with more than 1000 reads
#76 samples in the RRMS group with more than 1000 reads
###Since there are no samples with less than 1000 reads, there is no need to further prune

#Extracting family data from phyloseq object
family = tax_glom(MS_phyloseq,'Family')
ntaxa(MS_phyloseq); ntaxa(family)
# At the OTU level, there are 2752 taxa
#At the family level, there are 99 taxa

#Extracting family data from phyloseq object
genus = tax_glom(MS_phyloseq,'Genus')
ntaxa(MS_phyloseq); ntaxa(genus)
# At the OTU level, there are 2752 taxa
#At the genus level, there are 299 taxa

## Additional filtering
# First you define a function designed to work on a vector. For each x in the input vector, x will be divided by the sum of all x's in the vector.
#Takes the ASV abundance in a specific family and divides it by the total within that sample
calculate_relative_abundance <- function(x) x / sum(x)

# We'll only include things that are at least 0.1% abundant (0.001) across all samples
#Filtering out any ASVs with an abundance lower than 0.01%
total_counts <- taxa_sums(family) # Total sum for that taxa
relative_abundance <- calculate_relative_abundance(total_counts) # overall proportion of each bug
abundant <- relative_abundance > 0.001 # is each bug above the threshold? TRUE if so.
family <- prune_taxa(abundant, family) # subsetting our phyloseq object - Take only bugs above threshold
family 
# Now there are only 41 taxa, compared to 99 before filtering, and there are 151 samples

###ALDEX2 Analysis
library(ALDEx2)
set.seed(421)
#This analysis is for disease_course (comparing MS to control)
s = family@sam_data %>% as.matrix() %>% as.data.frame()
m = model.matrix(~ disease_course, data = s)
o = family@otu_table %>% as.matrix() %>% as.data.frame()

x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x)

colnames(df)

# Filter for significant features
significant_features <- df %>%
  filter(`disease_courseRRMS:pval` < 0.05) %>%
  arrange(desc(abs(`disease_courseRRMS:Est`))) 

# View the top results
head(significant_features)

# Add row names (ASV identifiers) as a new column
df_with_asv <- significant_features %>%
  rownames_to_column(var = "ASV")

# View the updated table with ASV column
head(df_with_asv)

# Extract taxonomy table from phyloseq object
taxonomy <- tax_table(MS_phyloseq)
taxonomy_df <- as.data.frame(taxonomy)

# Add row names as a new column in the taxonomy_df
taxonomy_df$`Feature_ID` <- rownames(taxonomy_df)
# View the updated taxonomy_df with the row names added as 'Feature ID'
head(taxonomy_df)

# Get the ASVs from the ALDEx2 results (assuming they are stored in df$ASV)
asv_ids <- df_with_asv$ASV

# Subset the taxonomy table to match the ASVs
taxonomy_sub <- taxonomy_df[taxonomy_df$`Feature_ID` %in% asv_ids, ]

# View the matched taxonomy information
head(taxonomy_sub)

# Merge the significant features with taxonomic information (e.g., Family)
significant_taxonomy <- merge(df_with_asv, taxonomy_sub, by.x = "ASV", by.y = "Feature_ID")

# View the merged data
head(significant_taxonomy)

library(ggplot2)

# Create the bar plot for significant features
MS_Control <- ggplot(significant_taxonomy, aes(x = `disease_courseRRMS:Est`, y = reorder(Genus, `disease_courseRRMS:Est`))) +
  geom_bar(stat = "identity") +  # Plot actual log2 fold change values
  geom_errorbar(aes(xmin = `disease_courseRRMS:Est` - `disease_courseRRMS:SE`, xmax = `disease_courseRRMS:Est` + `disease_courseRRMS:SE`), width = 0.2) +  # Add error bars
  theme_minimal() +  # Clean theme
  coord_flip() +  # Flip coordinates to make it easier to read
  labs(x = "log2 Fold Change", y = "Genus", title = "MS vs Control: Significant log2 Fold Changes by Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_blank())  # Rotate x-axis labels for readability

# ggsave saves the last-generated plot
ggsave(filename="MS_control.png",MS_Control)

###ALDEX2 Analysis part 2
library(ALDEx2)
set.seed(421)
#This analysis is for MS_asthma (comparing MS+asthma to control)
s1 = family@sam_data %>% as.matrix() %>% as.data.frame()
m1 = model.matrix(~ MS_asthma, data = s1)
o1 = family@otu_table %>% as.matrix() %>% as.data.frame()

x1 = aldex.clr(o1,m1,mc.samples=128)
df1 = aldex.glm(x1)

colnames(df1)

# Filter for significant features
significant_features1 <- df1 %>%
  filter(`MS_asthmaYes:pval` < 0.05) %>%
  arrange(desc(abs(`MS_asthmaYes:Est`))) 

# View the top results
head(significant_features1)

# Add row names (ASV identifiers) as a new column
df_with_asv1 <- significant_features1 %>%
  rownames_to_column(var = "ASV")

# View the updated table with ASV column
head(df_with_asv1)

# Extract taxonomy table from phyloseq object
taxonomy <- tax_table(MS_phyloseq)
taxonomy_df <- as.data.frame(taxonomy)
# Add row names as a new column in the taxonomy_df
taxonomy_df$`Feature_ID` <- rownames(taxonomy_df)
# View the updated taxonomy_df with the row names added as 'Feature ID'
head(taxonomy_df)

# Get the ASVs from the ALDEx2 results (assuming they are stored in df$ASV)
asv_ids1 <- df_with_asv1$ASV

# Subset the taxonomy table to match the ASVs
taxonomy_sub1 <- taxonomy_df[taxonomy_df$`Feature_ID` %in% asv_ids1, ]

# View the matched taxonomy information
head(taxonomy_sub1)

# Merge the significant features with taxonomic information (e.g., Family)
significant_taxonomy1 <- merge(df_with_asv1, taxonomy_sub1, by.x = "ASV", by.y = "Feature_ID")

# View the merged data
head(significant_taxonomy1)

library(ggplot2)

# Create the bar plot for significant features
MS_Asthma_vs_Control <- ggplot(significant_taxonomy1, aes(x = `MS_asthmaYes:Est`, y = reorder(Genus, `MS_asthmaYes:Est`))) +
  geom_bar(stat = "identity") +  # Plot actual log2 fold change values
  geom_errorbar(aes(xmin = `MS_asthmaYes:Est` - `MS_asthmaYes:SE`, xmax = `MS_asthmaYes:Est` + `MS_asthmaYes:SE`), width = 0.2) +  # Add error bars
  theme_minimal() +  # Clean theme
  coord_flip() +  # Flip coordinates to make it easier to read
  labs(x = "log2 Fold Change", y = "Genus", title = "MS_Asthma vs Control: Significant log2 Fold Changes by Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_blank())  # Rotate x-axis labels for readability

# ggsave saves the last-generated plot
ggsave(filename="MS_Asthma vs Control.png",MS_Asthma_vs_Control)

###ALDEX2 Analysis part 3
library(ALDEx2)
set.seed(421)
#This analysis is for asthma vs control
s2 = family@sam_data %>% as.matrix() %>% as.data.frame()
m2 = model.matrix(~ asthma, data = s2)
o2 = family@otu_table %>% as.matrix() %>% as.data.frame()

x2 = aldex.clr(o2,m2,mc.samples=128)
df2 = aldex.glm(x2)

colnames(df2)

# Filter for significant features
significant_features2 <- df2 %>%
  filter(`asthma1:pval` < 0.05) %>%
  arrange(desc(abs(`asthma1:Est`))) 

# View the top results
head(significant_features2)

# Add row names (ASV identifiers) as a new column
df_with_asv2 <- significant_features2 %>%
  rownames_to_column(var = "ASV")

# View the updated table with ASV column
head(df_with_asv2)

# Extract taxonomy table from phyloseq object
taxonomy <- tax_table(MS_phyloseq)
taxonomy_df <- as.data.frame(taxonomy)
# Add row names as a new column in the taxonomy_df
taxonomy_df$`Feature_ID` <- rownames(taxonomy_df)
# View the updated taxonomy_df with the row names added as 'Feature ID'
head(taxonomy_df)

# Get the ASVs from the ALDEx2 results (assuming they are stored in df$ASV)
asv_ids2 <- df_with_asv2$ASV

# Subset the taxonomy table to match the ASVs
taxonomy_sub2 <- taxonomy_df[taxonomy_df$`Feature_ID` %in% asv_ids2, ]

# View the matched taxonomy information
head(taxonomy_sub2)

# Merge the significant features with taxonomic information (e.g., Family)
significant_taxonomy2 <- merge(df_with_asv2, taxonomy_sub2, by.x = "ASV", by.y = "Feature_ID")

# View the merged data
head(significant_taxonomy2)

library(ggplot2)

# Create the bar plot for significant features
Asthma_vs_Control <- ggplot(significant_taxonomy2, aes(x = `asthma1:Est`, y = reorder(Genus, `asthma1:Est`))) +
  geom_bar(stat = "identity") +  # Plot actual log2 fold change values
  geom_errorbar(aes(xmin = `asthma1:Est` - `asthma1:SE`, xmax = `asthma1:Est` + `asthma1:SE`), width = 0.2) +  # Add error bars
  theme_minimal() +  # Clean theme
  coord_flip() +  # Flip coordinates to make it easier to read
  labs(x = "log2 Fold Change", y = "Genus", title = "Asthma vs Control: Significant log2 Fold Changes by Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_blank())  # Rotate x-axis labels for readability

# ggsave saves the last-generated plot
ggsave(filename="Asthma vs Control.png",Asthma_vs_Control)


