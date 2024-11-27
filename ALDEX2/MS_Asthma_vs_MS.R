#loading in required libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(DESeq2)


###ALDEX2 Analysis part 5
#Need to generate a new phyloseq object with altered metadata
#load in objects required for phyloseq object
metafp <- "ms_metadata_final_2.tsv"
metadata <- read_delim(metafp, delim="\t")

metadata$MS_asthma_vs_asthma <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 1, "MS_asthma", 
                                       ifelse(metadata$asthma == 1, "asthma", "No"))

metadata$MS_asthma_vs_MS <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 1, "MS_asthma", 
                                   ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 0, "MS", NA))

metadata$MS_vs_asthma <- ifelse(metadata$disease_course == "RRMS" & metadata$asthma == 0, "MS", 
                                ifelse(metadata$disease_course == "Control" & metadata$asthma == 1, "asthma", NA))

# Step 2: Filter the metadata to exclude rows where MS_asthma_vs_asthma is "No"
filtered_metadata <- metadata %>%
  filter(MS_asthma_vs_MS != "No")
filtered_metadata

otuFP <- "phyloseq/feature-table.txt"
otu <- read.delim(file=otuFP, skip=1)
otu

print(colnames(otu))

# Transform column names
colnames(otu) <- gsub("X(\\d+)\\.(\\d+)", "\\1-\\2", colnames(otu))
print(colnames(otu))

taxonomy <- read.delim('phyloseq/taxonomy.tsv')
taxonomy

phylo_tree <- read.tree('phyloseq/tree.nwk')
phylo_tree
class(phylo_tree)

##formatting otu matrix
#saving everything except for the first line into a matrix
otu_mat <- as.matrix(otu[,-1])
otu_mat

#Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`X.OTU.ID`
otu_mat


# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)
OTU

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(filtered_metadata[,-1])
samp_df 
# Make sampleids the rownames
rownames(samp_df)<- filtered_metadata$'sample-id'

# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)
SAMP

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- taxonomy %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
tax_mat

# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]

# Make sampleids the rownames
rownames(tax_mat) <- taxonomy$`Feature.ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)
TAX

#### Creating a phyloseq object ####
MS_phyloseq2 <- phyloseq(OTU, SAMP, TAX, phylo_tree)

# Let's see how our sequencing depth looks
hist(sample_sums(MS_phyloseq2),breaks = 30) 
table(below_1000 = sample_sums(MS_phyloseq2)<=1000) 
#  expt = MS_phyloseq@sam_data$disease_course) # shows us how many samples are below 1000
#All 151 samples have more than 1000 reads --> no need to further prune
#75 samples in the Control group (no MS) with more than 1000 reads
#76 samples in the RRMS group with more than 1000 reads
###Since there are no samples with less than 1000 reads, there is no need to further prune

#Extracting family data from phyloseq object
family = tax_glom(MS_phyloseq2,'Family')
ntaxa(MS_phyloseq); ntaxa(family)
# At the OTU level, there are 2752 taxa
#At the family level, there are 99 taxa

#Extracting family data from phyloseq object
genus = tax_glom(MS_phyloseq2,'Genus')
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
# Now there are only 40 taxa, compared to 99 before filtering, and there are 75 samples


library(ALDEx2)
set.seed(421)
#This analysis is for MS_asthma vs asthma
s = family@sam_data %>% as.matrix() %>% as.data.frame()
m = model.matrix(~ MS_asthma_vs_MS, data = s)
o = family@otu_table %>% as.matrix() %>% as.data.frame()

x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x)

colnames(df)

# Filter for significant features
significant_features <- df %>%
  filter(`MS_asthma_vs_MSMS_asthma:pval` < 0.05) %>%
  arrange(desc(abs(`MS_asthma_vs_MSMS_asthma:Est`))) 

# View the top results
head(significant_features)

# Add row names (ASV identifiers) as a new column
df_with_asv <- significant_features %>%
  rownames_to_column(var = "ASV")

# View the updated table with ASV column
head(df_with_asv)

# Extract taxonomy table from phyloseq object
taxonomy <- tax_table(MS_phyloseq2)
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
MS_asthma_vs_asthma <- ggplot(significant_taxonomy, aes(x = `MS_asthma_vs_asthmaMS_asthma:Est`, y = reorder(Genus, `MS_asthma_vs_asthmaMS_asthma:Est`))) +
  geom_bar(stat = "identity") +  # Plot actual log2 fold change values
  geom_errorbar(aes(xmin = `MS_asthma_vs_asthmaMS_asthma:Est` - `MS_asthma_vs_asthmaMS_asthma:SE`, xmax = `MS_asthma_vs_asthmaMS_asthma:Est` + `MS_asthma_vs_asthmaMS_asthma:SE`), width = 0.2) +  # Add error bars
  theme_minimal() +  # Clean theme
  coord_flip() +  # Flip coordinates to make it easier to read
  labs(x = "log2 Fold Change", y = "Genus", title = "MS_Asthma vs Asthma: Significant log2 Fold Changes by Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_blank())  # Rotate x-axis labels for readability

# ggsave saves the last-generated plot
ggsave(filename="MS Asthma vs Asthma.png",MS_asthma_vs_asthma)