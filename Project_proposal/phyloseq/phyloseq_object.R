#loading in required libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)

#load in objects required for phyloseq object
metafp <- "phyloseq/ms_metadata_final.tsv"
metadata <- read_delim(metafp, delim="\t")
metadata

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
samp_df <- as.data.frame(metadata[,-1])
samp_df 
# Make sampleids the rownames
rownames(samp_df)<- metadata$'sample-id'

# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)
SAMP

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- taxonomy_data %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
tax_mat

# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
tax_mat

# Make sampleids the rownames
rownames(tax_mat) <- taxonomy_data$`Sample Data`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)
TAX

#### Create phyloseq object ####
# Merge all into a phyloseq object
mouse_phyloseq <- phyloseq(OTU, SAMP, TAX, phylo_tree)

OTU 
