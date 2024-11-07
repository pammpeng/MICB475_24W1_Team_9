library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Beta diversity #####
bc_dm <- distance(ms_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(mpt_rare, method="PCoA", distance=bc_dm)

plot_ordination(mpt_rare, pcoa_bc, color = "body.site", shape="subject")

gg_pcoa <- plot_ordination(mpt_rare, pcoa_bc, color = "body.site", shape="subject") +
  labs(pch="Subject #", col = "Body Site")
gg_pcoa

ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)