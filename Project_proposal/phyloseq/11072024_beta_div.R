# Load the necessary package
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Beta diversity #####
bc_dm <- distance(ms_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(ms_rare, method="PCoA", distance=bc_dm)

plot_ordination(ms_rare, pcoa_bc, color = "disease_course", shape="asthma")

gg_pcoa <- plot_ordination(ms_rare, pcoa_bc, color = "disease_course", shape="asthma") +
  labs(pch="MS Status", col = "Asthma Status")
gg_pcoa


ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)

# Extract the distance matrix (e.g., Bray-Curtis)
dist_matrix <- distance(ms_rare, method = "bray")

# Extract the grouping variable from sample_data
group <- sample_data(ms_rare)$group  # Replace 'group' with your actual grouping variable name

# Perform PerMANOVA using adonis2 (since adonis is deprecated)
permanova_result <- adonis2(dist_matrix ~ group)

# View the results
print(permanova_result)
