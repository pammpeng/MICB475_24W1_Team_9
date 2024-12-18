library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Loading phyloseq object
load("MS_phyloseq.RData")

#### DESeq ####
#Generating a deseq object
MS_deseq <- phyloseq_to_deseq2(MS_phyloseq, ~`Condition`)
DESEQ_MS <- DESeq(MS_deseq)

## Since a zeros error was received in the previous lines of code, a '1' count is being added to all reads, then creating the deseq object.
MS_plus1 <- transform_sample_counts(MS_phyloseq, function(x) x+1)
MS_deseq <- phyloseq_to_deseq2(MS_plus1, ~`Condition`)
DESEQ_MS <- DESeq(MS_deseq)

###Generating results for the MS vs Control comparison by contrasting the two conditions, ensuring that "Control" is the reference group
res_MS_ctrl <- results(DESEQ_MS, tidy=TRUE, 
                       contrast = c("Condition","MS","Control"))
## Generating a table of results that associates the ASVs with their taxonomic information, log2foldchange, and any significant changes.
## Filtering for a p-adjusted value of less than 0.01 and a log2foldchange greater than 2
sigASVs_MS_ctrl <- res_MS_ctrl %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MS_ctrl)
## Returning only asv names
sigASVs_vec_MS_ctrl <- sigASVs_MS_ctrl %>%
  pull(ASV)

## Pruning phyloseq file to sort significant ASVs at the genus level
MS_DESeq <- prune_taxa(sigASVs_vec_MS_ctrl,MS_phyloseq)
sigASVs_MS_ctrl <- tax_table(MS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MS_ctrl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

## Generating a bar plot representing the genera that have a significant (p-adj<0.01) log2fold change
bar_plot_res_MS_ctrl <- ggplot(sigASVs_MS_ctrl) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS vs Control") 
ggsave(filename="bar_plot_MS_ctrl.png",bar_plot_res_MS_ctrl)

###Generating results for the Asthma vs Control comparison by contrasting the two conditions, ensuring that "Control" is the reference group
res_asthma_ctrl <- results(DESEQ_MS, tidy=TRUE, 
                           contrast = c("Condition","Asthma","Control"))
## Generating a table of results that associates the ASVs with their taxonomic information, log2foldchange, and any significant changes.
## Filtering for a p-adjusted value of less than 0.01 and a log2foldchange greater than 2
sigASVs_asthma_ctrl <- res_asthma_ctrl %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_asthma_ctrl)
# Returning only asv names
sigASVs_vec_asthma_ctrl <- sigASVs_asthma_ctrl %>%
  pull(ASV)

## Pruning phyloseq file to sort significant ASVs at the genus level
MS_DESeq1 <- prune_taxa(sigASVs_vec_asthma_ctrl,MS_phyloseq)
sigASVs_asthma_ctrl <- tax_table(MS_DESeq1) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_asthma_ctrl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()
## Generating a bar plot representing the genera that have a significant (p-adj<0.01) log2fold change
bar_plot_res_asthma_ctrl <- ggplot(sigASVs_asthma_ctrl) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Asthma vs Control") 
ggsave(filename="bar_plot_asthma_ctrl.png",bar_plot_res_asthma_ctrl)

###Generating results for the Asthma + MS vs Control comparison by contrasting the two conditions, ensuring that "Control" is the reference group
res_MSasthma_ctrl <- results(DESEQ_MS, tidy=TRUE, 
                             contrast = c("Condition","Asthma_MS","Control"))
## Generating a table of results that associates the ASVs with their taxonomic information, log2foldchange, and any significant changes.
## Filtering for a p-adjusted value of less than 0.01 and a log2foldchange greater than 2
sigASVs_MSasthma_ctrl <- res_MSasthma_ctrl %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MSasthma_ctrl)
## Returning only asv names
sigASVs_vec_MSasthma_ctrl <- sigASVs_MSasthma_ctrl %>%
  pull(ASV)

## Pruning phyloseq file to sort significant ASVs at the genus level
MS_DESeq2 <- prune_taxa(sigASVs_vec_MSasthma_ctrl,MS_phyloseq)
sigASVs_MSasthma_ctrl <- tax_table(MS_DESeq2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MSasthma_ctrl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

## Generating a bar plot representing the genera that have a significant (p-adj<0.01) log2fold change
bar_plot_res_MSasthma_ctrl <- ggplot(sigASVs_MSasthma_ctrl) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS + Asthma vs Control") 
ggsave(filename="bar_plot_MSasthma_ctrl.png",bar_plot_res_MSasthma_ctrl)

###Generating results for the Asthma + MS vs MS comparison by contrasting the two conditions, ensuring that "MS" is the reference group
res_MSasthma_MS <- results(DESEQ_MS, tidy=TRUE, 
                           contrast = c("Condition","Asthma_MS","MS"))
## Generating a table of results that associates the ASVs with their taxonomic information, log2foldchange, and any significant changes.
## Filtering for a p-adjusted value of less than 0.01 and a log2foldchange greater than 2
sigASVs_MSasthma_MS <- res_MSasthma_MS %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MSasthma_MS)
# Returning only asv names
sigASVs_vec_MSasthma_MS <- sigASVs_MSasthma_MS %>%
  pull(ASV)

## Pruning phyloseq file to sort significant ASVs at the genus level
MS_DESeq3 <- prune_taxa(sigASVs_vec_MSasthma_MS,MS_phyloseq)
sigASVs_MSasthma_MS <- tax_table(MS_DESeq3) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MSasthma_MS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

## Generating a bar plot representing the genera that have a significant (p-adj<0.01) log2fold change
bar_plot_res_MSasthma_MS <- ggplot(sigASVs_MSasthma_MS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS + Asthma vs MS") 
ggsave(filename="bar_plot_MSasthma_MS.png",bar_plot_res_MSasthma_MS)

###Generating results for the Asthma + MS vs Asthma comparison by contrasting the two conditions, ensuring that "Asthma" is the reference group
res_MSasthma_asthma <- results(DESEQ_MS, tidy=TRUE, 
                               contrast = c("Condition","Asthma_MS","Asthma"))

## Generating a table of results that associates the ASVs with their taxonomic information, log2foldchange, and any significant changes.
## Filtering for a p-adjusted value of less than 0.01 and a log2foldchange greater than 2
sigASVs_MSasthma_asthma <- res_MSasthma_asthma %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MSasthma_asthma)
# Returning only asv names
sigASVs_vec_MSasthma_asthma <- sigASVs_MSasthma_asthma %>%
  pull(ASV)

## Pruning phyloseq file to sort significant ASVs at the genus level
MS_DESeq4 <- prune_taxa(sigASVs_vec_MSasthma_asthma,MS_phyloseq)
sigASVs_MSasthma_asthma <- tax_table(MS_DESeq4) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MSasthma_asthma) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

## Generating a bar plot representing the genera that have a significant (p-adj<0.01) log2fold change
bar_plot_res_MSasthma_asthma <- ggplot(sigASVs_MSasthma_asthma) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS + Asthma vs Asthma") 
ggsave(filename="bar_plot_MSasthma_asthma.png",bar_plot_res_MSasthma_asthma)

###Generating results for the MS vs Asthma comparison by contrasting the two conditions, ensuring that "Asthma" is the reference group
res_MS_asthma <- results(DESEQ_MS, tidy=TRUE, 
                         contrast = c("Condition","MS","Asthma"))
## Generating a table of results that associates the ASVs with their taxonomic information, log2foldchange, and any significant changes.
## Filtering for a p-adjusted value of less than 0.01 and a log2foldchange greater than 2
sigASVs_MS_asthma <- res_MS_asthma %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MS_asthma)
# Returning only asv names
sigASVs_vec_MS_asthma <- sigASVs_MS_asthma %>%
  pull(ASV)

## Pruning phyloseq file to sort significant ASVs at the genus level
MS_DESeq5 <- prune_taxa(sigASVs_vec_MS_asthma,MS_phyloseq)
sigASVs_MS_asthma <- tax_table(MS_DESeq5) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MS_asthma) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

## Generating a bar plot representing the genera that have a significant (p-adj<0.01) log2fold change
bar_plot_res_MS_asthma <- ggplot(sigASVs_MS_asthma) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS vs Asthma") 
ggsave(filename="bar_plot_MS_asthma.png",bar_plot_res_MS_asthma)

### Note that the graphs generated above are not included in the manuscript or the GitHub as they were used for a preliminary analysis of the data.
### Further filtering of the data was applied to select graphs as detailed below to find genera that were present in all comparisons (except MS+Asthma vs MS).

### Determining genera that are present in all five comparisons of interest (everything except MS+Asthma vs Asthma was included as the 
### number of  genera that increased and decreased in abundance was similar).

## Extracting genera that significantly change in log2FoldChange (p-adj<0.01)
extract_genera <- function(DESEQ_MS, contrast) {
  # Generating the results table for each specified contrast
  res <- results(DESEQ_MS, tidy = TRUE, contrast = contrast)
 
  ## Filtering for significantly differentially abundant ASVs (adjusted p-value < 0.01 and |log2FoldChange| > 2)
  sigASVs <- res %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
    dplyr::rename(ASV = row) %>%
    pull(ASV)
 
   ## Pruning phyloseq file to sort significant ASVs at a taxonomic level
  pruned_phyloseq <- prune_taxa(sigASVs, MS_phyloseq)
  taxa_info <- tax_table(pruned_phyloseq) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ASV")
  
  # Returning unique genera associated with the ASVs
  unique(taxa_info$Genus)
}

# Extracting the significantly changed genera for each comparison of interest
genera_MS_Control <- extract_genera(DESEQ_MS, c("Condition", "MS", "Control"))
genera_Asthma_Control <- extract_genera(DESEQ_MS, c("Condition", "Asthma", "Control"))
genera_MS_Asthma <- extract_genera(DESEQ_MS, c("Condition", "MS", "Asthma"))
genera_MSasthma_Control <- extract_genera(DESEQ_MS, c("Condition", "Asthma_MS", "Control"))
genera_MSasthma_Asthma <- extract_genera(DESEQ_MS, c("Condition", "Asthma_MS", "Asthma"))


# Identifying the common genera present in all comparisons
common_genera <- Reduce(intersect, list(
  genera_MS_Control,
  genera_Asthma_Control,
  genera_MS_Asthma,
  genera_MSasthma_Control,
  genera_MSasthma_Asthma
)) %>% setdiff("g__Bacteroides")
#Note that Bacteroides is excluded because it is not actually present in the MS vs Asthma comparison.

common_genera

# Function to process results and generate bar plots
process_results <- function(DESEQ_MS, contrast, title, prefix) {
  # Obtain results for each specified contrast
  res <- results(DESEQ_MS, tidy = TRUE, contrast = contrast)
  
  # Filtering significant differentially abundant ASVs
  sigASVs <- res %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
    dplyr::rename(ASV = row)
  
  # Subsetting the phyloseq object to include only significant ASVs and joining with taxonomy data
  pruned_phyloseq <- prune_taxa(sigASVs$ASV, MS_phyloseq)
  taxa_info <- tax_table(pruned_phyloseq) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ASV")
  
  ## Merging ASV data with taxonomic information and formatting for plotting
  sigASVs <- sigASVs %>%
    left_join(taxa_info, by = "ASV") %>%
    arrange(log2FoldChange) %>%
    mutate(Genus = make.unique(Genus)) %>%
    mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
    drop_na()
  
  # Labelling genera as either "common" (present in all comparisons) or "not common"
  sigASVs <- sigASVs %>%
    mutate(`Common genera` = ifelse(Genus %in% common_genera, "Yes", "No"))
  
  # Generating a bar plot showing log2FoldChange for each genus, with common genera in red, and others in grey
  bar_plot <- ggplot(sigASVs) +
    geom_bar(aes(x = Genus, y = log2FoldChange, fill = `Common genera`), stat = "identity") +
    geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle(title) + theme(
      plot.background = element_rect(
        fill = "white",
        colour = "white"))
  
  # Saving the bar plot
  ggsave(filename = paste0(prefix, "_bar_plot.png"), bar_plot, width = 8, height = 6)
  
  # Returning the processed results
  return(sigASVs)
}

# Defining all pairwise comparisons of interest so that plots can be generated for each.
comparisons <- list(
  list(contrast = c("Condition", "MS", "Control"), title = "MS vs Control", prefix = "MS_Control"),
  list(contrast = c("Condition", "Asthma", "Control"), title = "Asthma vs Control", prefix = "Asthma_Control"),
  list(contrast = c("Condition", "MS", "Asthma"), title = "MS vs Asthma", prefix = "MS_Asthma"),
  list(contrast = c("Condition", "Asthma_MS", "Control"), title = "MS+Asthma vs Control", prefix = "MSasthma_Control"),
  list(contrast = c("Condition", "Asthma_MS", "Asthma"), title = "MS+Asthma vs Asthma", prefix = "MSasthma_Asthma")
)

# Applying the updated process_results function to each comparison
results_list_final <- lapply(comparisons, function(comp) {
  process_results(DESEQ_MS, comp$contrast, comp$title, comp$prefix)
})


