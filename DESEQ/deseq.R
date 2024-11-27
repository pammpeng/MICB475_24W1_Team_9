library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("MS_phyloseq.RData")

#### DESeq ####
## NOTE: If you get a zeros error, then you need to add '1' count to all reads
MS_plus1 <- transform_sample_counts(MS_phyloseq, function(x) x+1)
MS_deseq <- phyloseq_to_deseq2(MS_plus1, ~`Condition`)
DESEQ_MS <- DESeq(MS_deseq)

###Results for MS vs Control
res_MS_ctrl <- results(DESEQ_MS, tidy=TRUE, 
               #this will ensure that Control is your reference group
               contrast = c("Condition","MS","Control"))
#Volcano plot for MS vs Control
vol_plot_res_MS_ctrl  <- res_MS_ctrl  %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("MS vs Control") 
vol_plot_res_MS_ctrl
ggsave(filename="vol_plot_MS_ctrl.png",vol_plot_res_MS_ctrl)

# To get table of results
sigASVs_MS_ctrl <- res_MS_ctrl %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MS_ctrl)
# Get only asv names
sigASVs_vec_MS_ctrl <- sigASVs_MS_ctrl %>%
  pull(ASV)

# Prune phyloseq file
MS_DESeq <- prune_taxa(sigASVs_vec_MS_ctrl,MS_phyloseq)
sigASVs_MS_ctrl <- tax_table(MS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MS_ctrl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

bar_plot_res_MS_ctrl <- ggplot(sigASVs_MS_ctrl) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS vs Control") 
ggsave(filename="bar_plot_MS_ctrl.png",bar_plot_res_MS_ctrl)

###Results for Asthma vs Control
res_asthma_ctrl <- results(DESEQ_MS, tidy=TRUE, 
               contrast = c("Condition","Asthma","Control"))
#Volcano plot for Asthma vs Control
vol_plot_res_asthma_ctrl  <- res_asthma_ctrl  %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Asthma vs Control")
vol_plot_res_asthma_ctrl
ggsave(filename="vol_plot_asthma_ctrl.png",vol_plot_res_asthma_ctrl)

# To get table of results
sigASVs_asthma_ctrl <- res_asthma_ctrl %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_asthma_ctrl)
# Get only asv names
sigASVs_vec_asthma_ctrl <- sigASVs_asthma_ctrl %>%
  pull(ASV)

# Prune phyloseq file
MS_DESeq1 <- prune_taxa(sigASVs_vec_asthma_ctrl,MS_phyloseq)
sigASVs_asthma_ctrl <- tax_table(MS_DESeq1) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_asthma_ctrl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

bar_plot_res_asthma_ctrl <- ggplot(sigASVs_asthma_ctrl) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Asthma vs Control") 
ggsave(filename="bar_plot_asthma_ctrl.png",bar_plot_res_asthma_ctrl)

###Results for MS+Asthma vs Control
res_MSasthma_ctrl <- results(DESEQ_MS, tidy=TRUE, 
                           contrast = c("Condition","Asthma_MS","Control"))
#Volcano plot for MS+Asthma vs Control
vol_plot_res_MSasthma_ctrl  <- res_MSasthma_ctrl  %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("MS+Asthma vs Control")
vol_plot_res_MSasthma_ctrl
ggsave(filename="vol_plot_MSasthma_ctrl.png",vol_plot_res_MSasthma_ctrl)

# To get table of results
sigASVs_MSasthma_ctrl <- res_MSasthma_ctrl %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MSasthma_ctrl)
# Get only asv names
sigASVs_vec_MSasthma_ctrl <- sigASVs_MSasthma_ctrl %>%
  pull(ASV)

# Prune phyloseq file
MS_DESeq2 <- prune_taxa(sigASVs_vec_MSasthma_ctrl,MS_phyloseq)
sigASVs_MSasthma_ctrl <- tax_table(MS_DESeq2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MSasthma_ctrl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

bar_plot_res_MSasthma_ctrl <- ggplot(sigASVs_MSasthma_ctrl) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS + Asthma vs Control") 
ggsave(filename="bar_plot_MSasthma_ctrl.png",bar_plot_res_MSasthma_ctrl)

###Results for MS+Asthma vs MS
res_MSasthma_MS <- results(DESEQ_MS, tidy=TRUE, 
                             contrast = c("Condition","Asthma_MS","MS"))
#Volcano plot for MS+Asthma vs MS
vol_plot_res_MSasthma_MS  <- res_MSasthma_MS  %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("MS+Asthma vs MS")
vol_plot_res_MSasthma_MS
ggsave(filename="vol_plot_MSasthma_MS.png",vol_plot_res_MSasthma_MS)

# To get table of results
sigASVs_MSasthma_MS <- res_MSasthma_MS %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MSasthma_MS)
# Get only asv names
sigASVs_vec_MSasthma_MS <- sigASVs_MSasthma_MS %>%
  pull(ASV)

# Prune phyloseq file
MS_DESeq3 <- prune_taxa(sigASVs_vec_MSasthma_MS,MS_phyloseq)
sigASVs_MSasthma_MS <- tax_table(MS_DESeq3) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MSasthma_MS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

bar_plot_res_MSasthma_MS <- ggplot(sigASVs_MSasthma_MS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS + Asthma vs MS") 
ggsave(filename="bar_plot_MSasthma_MS.png",bar_plot_res_MSasthma_MS)

###Results for MS+Asthma vs Asthma
res_MSasthma_asthma <- results(DESEQ_MS, tidy=TRUE, 
                             contrast = c("Condition","Asthma_MS","Asthma"))
#Volcano plot for MS+Asthma vs Asthma
vol_plot_res_MSasthma_asthma  <- res_MSasthma_asthma  %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("MS+Asthma vs Asthma")
vol_plot_res_MSasthma_asthma
ggsave(filename="vol_plot_MSasthma_Asthma.png",vol_plot_res_MSasthma_asthma)

# To get table of results
sigASVs_MSasthma_asthma <- res_MSasthma_asthma %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MSasthma_asthma)
# Get only asv names
sigASVs_vec_MSasthma_asthma <- sigASVs_MSasthma_asthma %>%
  pull(ASV)

# Prune phyloseq file
MS_DESeq4 <- prune_taxa(sigASVs_vec_MSasthma_asthma,MS_phyloseq)
sigASVs_MSasthma_asthma <- tax_table(MS_DESeq4) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MSasthma_asthma) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

bar_plot_res_MSasthma_asthma <- ggplot(sigASVs_MSasthma_asthma) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS + Asthma vs Asthma") 
ggsave(filename="bar_plot_MSasthma_asthma.png",bar_plot_res_MSasthma_asthma)

###Results for MS vs Asthma
res_MS_asthma <- results(DESEQ_MS, tidy=TRUE, 
                             contrast = c("Condition","MS","Asthma"))
#Volcano plot for MS vs Asthma
vol_plot_res_MS_asthma  <- res_MS_asthma  %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("MS vs Asthma")
vol_plot_res_MS_asthma
ggsave(filename="vol_plot_MS_Asthma.png",vol_plot_res_MS_asthma)

# To get table of results
sigASVs_MS_asthma <- res_MS_asthma %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_MS_asthma)
# Get only asv names
sigASVs_vec_MS_asthma <- sigASVs_MS_asthma %>%
  pull(ASV)

# Prune phyloseq file
MS_DESeq5 <- prune_taxa(sigASVs_vec_MS_asthma,MS_phyloseq)
sigASVs_MS_asthma <- tax_table(MS_DESeq5) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_MS_asthma) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus))) %>%
  drop_na()

bar_plot_res_MS_asthma <- ggplot(sigASVs_MS_asthma) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("MS vs Asthma") 
ggsave(filename="bar_plot_MS_asthma.png",bar_plot_res_MS_asthma)









