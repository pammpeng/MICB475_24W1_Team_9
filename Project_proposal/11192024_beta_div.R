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

head(sample_data(ms_rare))

plot_ordination(ms_rare, pcoa_bc, color = "group")

gg_pcoa <- plot_ordination(ms_rare, pcoa_bc, color = "disease_course", shape="asthma") +
  labs(pch="Asthma Status", col = "MS Status")+
  stat_ellipse(aes(group = group), level = 0.95) 
gg_pcoa


ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)

# Extract the distance matrix (e.g., Bray-Curtis)
dist_matrix <- distance(ms_rare, method = "bray")

# Extract the grouping variable from sample_data
asthma <- sample_data(ms_rare)$asthma  # Replace 'group' with your actual grouping variable name

# Perform PerMANOVA using adonis2 (since adonis is deprecated)
permanova_result <- adonis2(dist_matrix ~ asthma)

# View the results
print(permanova_result)

alphadiv <- estimate_richness(ms_rare)
samp_dat <- sample_data(ms_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

ml_ms_asthma <- lm(Shannon ~ `disease_course`*`asthma`, data=samp_dat_wdiv)
summary(aov(ml_ms_asthma))
TukeyHSD(aov(ml_ms_asthma))

samp_data_wdiv_plot <- ggplot(samp_dat_wdiv) + geom_boxplot(aes(x=disease_course, y=Shannon)) +
  facet_grid(~factor(`asthma`))+
  theme_classic()
summary(ml_ms_asthma)
ggsave("samp_data_div_plot.png", width = 8, height = 6, dpi = 300)

head(sample_data(ms_rare))
sample_data(ms_rare) <- sample_data(ms_rare)[, -8]

dm_braycurtis <- vegdist(t(otu_table(ms_rare)), method="bray")  # Bray-Curtis distance
dat <- as.data.frame(sample_data(ms_rare))
dm_braycurtis

dat$disease_course <- as.factor(dat$disease_course)
class(dat$asthma)
class(dat$disease_course)

adonis2(dm_unifrac ~ disease_course*asthma, data=dat)                         


class(dat)


# Ensure 'dat' is a data frame
class(dat)  # Should return "data.frame"
dat <- data.frame(sample_data(ms_rare))  # If it's not a data frame, convert it
class(dat)

# make 'disease_course' and 'asthma' into factors
dat$disease_course <- as.factor(dat$disease_course)
dat$asthma <- as.factor(dat$asthma)


# If that works, try with the interaction term
adonis2(dm_braycurtis ~ disease_course * asthma, data = dat)
