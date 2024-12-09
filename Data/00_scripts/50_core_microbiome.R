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


# Create lists of the variables to compare 
condition_list_all <- list("MS + Asthma" = asthma_MS_ASVs, MS = noasthma_MS_ASVs, 
                            Healthy = control_ASVs, Asthma = asthma_noMS_ASVs)
condition_list_MSvsControl <- list(Healthy = control_ASVs,MS = noasthma_MS_ASVs)
condition_list_MSAsthmavsControl <- list(Healthy = control_ASVs,"MS + Asthma" = asthma_MS_ASVs)
condition_list_minus_asthma <- list(Healthy = control_ASVs,"MS + Asthma" = asthma_MS_ASVs, MS = noasthma_MS_ASVs)
condition_list_MS_asthmavsno <- list("MS + Asthma" = asthma_MS_ASVs, MS = noasthma_MS_ASVs)
condition_list_controlvsasthma <- list(Healthy = control_ASVs, Asthma = asthma_noMS_ASVs)
condition_list_asthmavsasmtha_MS <- list(Asthma = asthma_noMS_ASVs, "MS + Asthma" = asthma_MS_ASVs)

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
venn_controlvsasthma <- ggVennDiagram(x = condition_list_controlvsasthma) + coord_flip() + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))
venn_asthmavsasthma_MS <- ggVennDiagram(x = condition_list_asthmavsasmtha_MS) + coord_flip() + theme(
  plot.background = element_rect(
    fill = "white",
    colour = "white"))

# Save Venn diagrams
ggsave("core_microbiome_venn.png", venn_all, height = 8, width = 12)
ggsave("venn_MS_vs_Control.png", venn_MSvsControl, height = 6, width = 8)
ggsave("venn_MS_Asthma_vs_Control.png", venn_MSAsthmavsControl, height = 6, width = 8)
ggsave("venn_minus_Asthma.png", venn_minus_Asthma)
ggsave("ven_MS.png", venn_MS_asthmavsno)
ggsave("ven_control_vs_asthma.png", venn_controlvsasthma)
