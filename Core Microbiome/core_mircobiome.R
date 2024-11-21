#Load Packages 
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load Phyloseq Object ####
load("MS_phyloseq.RData")

#### Core Microbiome Analysis ####

# Convert absolute abundance to relative abundance
MS_RA <- transform_sample_counts(MS_phyloseq, fun=function(x) x/sum(x))

# Filter dataset by condition 
asthma_MS <- subset_samples(MS_RA, `Condition`=="Asthma_MS")
noasthma_MS <- subset_samples(MS_RA, `Condition`=="MS")
control <- subset_samples(MS_RA, `Condition`=="Control")
asthma_noMS <- subset_samples(MS_RA, `Condition`=="Asthma")

# Determining ASVs present in 50% of samples 
asthma_MS_ASVs <- core_members(asthma_MS, detection=0.01, prevalence = 0.5)
noasthma_MS_ASVs <- core_members(noasthma_MS, detection=0.001, prevalence = 0.5)
control_ASVs <- core_members(control, detection=0.001, prevalence = 0.5)
asthma_noMS_ASVs <- core_members(asthma_noMS, detection=0.001, prevalence = 0.5)

# Create plot of ASVs relative abundance 
condition_taxa_plot <- prune_taxa(noasthma_MS_ASVs,MS_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Condition`, scales ="free")

# Create lists of the variables to compare 
condition_list_all <- list(Asthma_MS = asthma_MS_ASVs, MS = noasthma_MS_ASVs, 
                            Control = control_ASVs, Asthma = asthma_noMS_ASVs)
condition_list_MSvsControl <- list(Control = control_ASVs,MS = noasthma_MS_ASVs)
condition_list_MSAsthmavsControl <- list(Control = control_ASVs,Asthma_MS = asthma_MS_ASVs)
condition_list_minus_asthma <- list(Control = control_ASVs,Asthma_MS = asthma_MS_ASVs, MS = noasthma_MS_ASVs)
condition_list_MS_asthmavsno <- list(Asthma_MS = asthma_MS_ASVs, MS = noasthma_MS_ASVs)

# Create a Venn diagrams 
venn_all <- ggVennDiagram(x = condition_list_all) + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))
venn_MSvsControl <- ggVennDiagram(x = condition_list_MSvsControl) + coord_flip() + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))
venn_MSAsthmavsControl <- ggVennDiagram(x = condition_list_MSAsthmavsControl) + coord_flip() + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))
venn_minus_Asthma <- ggVennDiagram(x = condition_list_minus_asthma) + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))
venn_MS_asthmavsno <- ggVennDiagram(x = condition_list_MS_asthmavsno) + coord_flip() + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))

# Save Venn diagrams
ggsave("venn_all.png", venn_all, height = 6, width = 10)
ggsave("venn_MS_vs_Control.png", venn_MSvsControl, height = 6, width = 8)
ggsave("venn_MS_Asthma_vs_Control.png", venn_MSAsthmavsControl, height = 6, width = 8)
ggsave("venn_minus_Asthma.png", venn_minus_Asthma)
ggsave("ven_MS.png", venn_MS_asthmavsno)
