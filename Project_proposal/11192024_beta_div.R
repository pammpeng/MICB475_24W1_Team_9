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

gg_pcoa <- plot_ordination(ms_rare, pcoa_bc, color = "group") +
  labs(pch="Asthma Status", col = "Health Status")
gg_pcoa


ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)

library(vegan)

# Perform PCA (or PCoA) on the distance matrix
pcoa_result <- cmdscale(bc_dm , k = 2)  # k = 2 for two principal coordinates

# Convert PCA results to a data frame
pcoa_df <- data.frame(pcoa_result)

# Add the categorical variables to the data frame
pcoa_df$disease_course <- factor(sample_data(ms_rare)$disease_course)
pcoa_df$asthma <- factor(sample_data(ms_rare)$asthma)

# Check the structure of the data
head(pcoa_df)

# Fit linear models for each principal coordinate
lm_PC1 <- lm(X1 ~ disease_course * asthma, data = pcoa_df)
lm_PC2 <- lm(X2 ~ disease_course * asthma, data = pcoa_df)

# Summary of the models
summary(lm_PC1)
summary(lm_PC2)

library(ggplot2)
pcoa_plot_2 <- ggplot(pcoa_df, aes(x = X1, y = X2, color = disease_course, shape = asthma)) +
  geom_point(size = 3) +
  labs(title = "PCoA Plot by MS and Asthma Status",
       x = "PC1", y = "PC2") +
  theme_classic()
pcoa_plot_2

ggsave("plot_pcoa_2.png"
       , pcoa_plot_2
       , height=4, width=5)

# Extract the distance matrix (e.g., Bray-Curtis)
dist_matrix <- distance(ms_rare, method = "bray")

# Extract the grouping variable from sample_data
asthma <- sample_data(ms_rare)$asthma  # Replace 'group' with your actual grouping variable name

# Perform PerMANOVA using adonis2 (since adonis is deprecated)
permanova_result <- adonis2(dist_matrix ~ asthma)


dm_braycurtis <- vegdist(t(otu_table(ms_rare)), method="bray")  # Bray-Curtis distance
dat <- as.data.frame(sample_data(ms_rare))
dm_braycurtis


dat$disease_course <- as.factor(dat$disease_course)
class(dat$asthma)
class(dat$disease_course)
class(get_sample(ms_rare))
ms_df <- as_data_frame(get_sample(ms_rare))
ms_df
adonis2(dm_braycurtis ~ disease_course*asthma, data=dat)                         


class(dat)


# Make dat into a dataframe
class(dat)  
dat <- data.frame(sample_data(ms_rare))
class(dat)

# make 'disease_course' and 'asthma' into factors
dat$disease_course <- as.factor(dat$disease_course)
dat$asthma <- as.factor(dat$asthma)


# Run Adonis2 to perform PERMANOVA
adonis2(dm_braycurtis ~ disease_course * asthma, data = dat)
