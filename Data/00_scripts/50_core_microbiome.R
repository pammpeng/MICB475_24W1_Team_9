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
asthma_MS_ASVs <- core_members(asthma_MS, detection=0.001, prevalence = 0.5)
noasthma_MS_ASVs <- core_members(noasthma_MS, detection=0.001, prevalence = 0.5)
control_ASVs <- core_members(control, detection=0.001, prevalence = 0.5)
asthma_noMS_ASVs <- core_members(asthma_noMS, detection=0.001, prevalence = 0.5)

# Create a list of all the variables to compare 
condition_list_all <- list("MS + Asthma" = asthma_MS_ASVs, MS = noasthma_MS_ASVs, 
                            Healthy = control_ASVs, Asthma = asthma_noMS_ASVs)

# Create a Venn diagram 
venn_all <- ggVennDiagram(x = condition_list_all) + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))

# Save Venn diagram
ggsave("core_microbiome_venn.png", venn_all, height = 8, width = 12)
