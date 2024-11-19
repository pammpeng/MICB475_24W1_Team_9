#loading in required libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(DESeq2)

#load in objects required for phyloseq object
metafp <- "phyloseq/ms_metadata_final_2.tsv"
metadata <- read_delim(metafp, delim="\t")

#Adding a new column to the metadata with the condition assignment (asthma, asthma+MS, etc)
metadata$Condition <- with(metadata, ifelse(disease_course == "RRMS" & asthma == 1, "Asthma_MS",
                                      ifelse(disease_course == "Control" & asthma == 1, "Asthma",
                                             ifelse(disease_course == "Control" & asthma == 0, "Control", 
                                                    ifelse(disease_course == "RRMS" & asthma == 0, "MS", NA)))))
metadata

#Making Condition a factor
metadata$Condition <- factor(metadata$Condition, levels = c("Control", "Asthma_MS", "Asthma", "MS"))

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
MS_phyloseq <- phyloseq(OTU, SAMP, TAX, phylo_tree)

#### DESeq ####

#Creating deseq
MS_deseq <- phyloseq_to_deseq2(MS_phyloseq, ~`Condition`)

## Adding '1' count to all reads
MS_plus1 <- transform_sample_counts(MS_phyloseq, function(x) x+1)
MS_deseq <- phyloseq_to_deseq2(MS_plus1, ~`Condition`)
DESEQ_MS <- DESeq(MS_deseq)

# Defining all contrasts
comparisons <- list(c("MS", "Control"),
                    c("Asthma", "Control"),
                    c("Asthma_MS", "Control"),
                    c("Asthma_MS", "MS"),
                    c("Asthma_MS", "Asthma"))

results_list <- list()

# Looping through each comparison
for (comp in comparisons) {
  res <- results(DESEQ_MS, contrast = c("Condition", comp[1], comp[2]))
  results_list[[paste(comp[1], "vs", comp[2], sep = "_")]] <- res
}
resultsNames(res)
View(res)

# Look at results 

#Creating a volcano plot
create_volcano_plot <- function(res, title = "Volcano Plot") {
  # Converting results to a data frame
  res_df <- as.data.frame(res)
  res_df$significant <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 2, "Significant", "Not Significant")
  
# Generating the plots
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(p-value)") +
    theme(legend.title = element_blank())
}

#Generating the plots in a loop
for (comp_name in names(results_list)) {
  res <- results_list[[comp_name]]
  
  # Create the volcano plot
  p <- create_volcano_plot(res, title = comp_name)
  
  # Print the plot
  print(p)
  
  # Save to file
  ggsave(filename = paste0("Volcano_Plot_", comp_name, ".png"), plot = p, width = 8, height = 6)
}

#Filtering for significant ASVs and generating plots
# List to store filtered significant ASVs and plots for each comparison
significant_ASVs <- list()
plots <- list()

# Loop through the results list
for (comp_name in names(results_list)) {
  # Extract DESeq2 results for the current comparison
  res <- results_list[[comp_name]]
  
  # Convert DESeq2 results to a data frame
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "ASV")
  
  # Filter for significant ASVs
  sigASVs_vec <- res_df %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
    pull(ASV)  # Get significant ASV names
  
  # Check if there are significant ASVs
  if (length(sigASVs_vec) == 0) {
    cat(paste0("No significant ASVs for comparison: ", comp_name, "\n"))
    next  # Skip to the next iteration
  }
  
  # Prune the phyloseq object for significant ASVs
  pruned_phyloseq <- prune_taxa(sigASVs_vec, MS_phyloseq)
  
  # Prepare taxonomy data
  taxonomy_df <- as.data.frame(as.matrix(tax_table(pruned_phyloseq))) %>%  # Ensure standard data frame
    rownames_to_column(var = "ASV")
  
  # Join taxonomy with filtered DESeq2 results
  sigASVs_df <- taxonomy_df %>%
    right_join(res_df, by = "ASV") %>%
    arrange(desc(abs(log2FoldChange))) %>%  # Sort by absolute log2 fold change
    head(20)  # Select only the top 20 ASVs
  
  # Store significant ASVs
  significant_ASVs[[comp_name]] <- sigASVs_df
  
  # Create the bar plot
  p <- ggplot(sigASVs_df) +
    geom_bar(aes(x = Genus, y = log2FoldChange, fill = Genus), stat = "identity") +
    geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Adjust text angle and size
      plot.margin = margin(1, 1, 2, 1, "cm")  # Add margins for better spacing
    ) +
    labs(
      title = paste("Top 20 Differentially Abundant Taxa for", comp_name),
      x = "Genus",
      y = "Log2 Fold Change"
    ) +
    scale_fill_discrete(name = "Genus")  # Add legend for genus groups
  
  # Save the plot in the list
  plots[[comp_name]] <- p
  
  # Print progress
  cat(paste0("Processed comparison: ", comp_name, "\n"))
}

# Save or view the plots
for (comp_name in names(plots)) {
  print(plots[[comp_name]])  # View in RStudio
  ggsave(filename = paste0("Barplot_", comp_name, ".png"), plot = plots[[comp_name]], width = 12, height = 10, dpi = 300)
}

# Summary: Print the number of significant ASVs for each comparison
for (comp_name in names(significant_ASVs)) {
  cat(paste0("Comparison: ", comp_name, " - ", nrow(significant_ASVs[[comp_name]]), " significant ASVs\n"))
}

