library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Beta diversity #####
bc_dm <- distance(ms_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(ms_rare, method="PCoA", distance=bc_dm)

plot_ordination(ms_rare, pcoa_bc, color = "asthma", shape = "disease_course")

gg_pcoa <- plot_ordination(ms_rare, pcoa_bc, color = "asthma", shape="disease_course") +
  labs(pch="Disease Status", col = "Asthma")+
  theme_minimal()
gg_pcoa

ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)
