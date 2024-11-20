library(tidyverse)
library(picante)

library(ggplot2)
library(dplyr)
library("ggpubr")

install.packages("ggpubr")


#### Alpha diversity ######
plot_richness(ms_rare)

# phylogenetic diversity


# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(ms_rare)), phy_tree(ms_rare),
                 include.root=F) 
?pd

# add PD to metadata table
sample_data(ms_rare)$PD <- phylo_dist$PD
sample_data(ms_rare)



sample_data(ms_rare)$group <- paste(sample_data(ms_rare)$disease_course, sample_data(ms_rare)$asthma, sep = "_")


# Extract the grouping variable from sample data
group_info <- sample_data(ms_rare)$group  # Replace 'group' with your actual grouping variable
group_info

# Create a data frame with Faith's PD and the group info
faith_pd_df <- data.frame(PD = sample_data(ms_rare)$PD, group = group_info)

# Check the first few rows of the data frame
head(faith_pd_df)
#Change group names
sample_data(ms_rare)$group <- gsub("Control_0", "Healthy", sample_data(ms_rare)$group)
sample_data(ms_rare)$group <- gsub("Control_1", "Asthma", sample_data(ms_rare)$group)
sample_data(ms_rare)$group <- gsub("RRMS_0", "MS", sample_data(ms_rare)$group)
sample_data(ms_rare)$group <- gsub("RRMS_1", "MS + Asthma", sample_data(ms_rare)$group)

# Extract the sample metadata from ms_rare
sample_metadata <- sample_data(ms_rare)

# Select the desired columns
new_faiths <- sample_metadata[, c("PD", "disease_course", "asthma")]

# View the new data frame
head(new_faiths)

sample_data(ms_rare)

head(new_faiths)
class(new_faiths)
new_faiths$PD <- as.numeric(new_faiths$PD)
new_faiths$disease_course <- factor(new_faiths$disease_course)
new_faiths$asthma <- factor(new_faiths$asthma)

#perform PERMANOVA test for multiple categorical variables
new_faiths
new_faiths <- data.frame(sample_data(new_faiths))
permanova_result <- adonis2(PD ~ disease_course * asthma, data = new_faiths, permutations = 999)
head(new_faiths)

colnames(new_faiths)

names(new_faiths)
permanova_result <- adonis2(as.formula("PD ~ disease_course * asthma"), data = new_faiths, permutations = 999)

# View the result
permanova_result

new_faiths <- data.frame(sample_data(new_faiths))
# Perform Kruskal-Wallis test to compare PD across groups
kruskal_test_result <- kruskal.test(PD ~ group, data = faith_pd_df)

# Print the Kruskal-Wallis test result
print(kruskal_test_result)

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ms_rare), aes(group, PD), fill = group) + 
  geom_boxplot() +
  xlab("Group") +
  ylab("Phylogenetic Diversity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# view plot
plot.pd

# Save the last plot as a PNG image
ggsave("plot_pd.png", width = 8, height = 6, dpi = 300)

#kruskal wallis test for pd 

# Extract group information from the sample data (e.g., disease status, treatment, etc.)
group_info <- sample_data(ms_rare)$group  # Adjust "group" to your actual variable
group_info

# Create a data frame combining Faith's PD and group information
faith_pd_df <- data.frame(faith_pd = faith_pd_values, group = group_info)

# View the first few rows of the data frame
head(faith_pd_df)
kruskal_test_result <- kruskal.test(PD ~ group, data = faith_pd_df)
kruskal_test_result


# Calculate alpha diversity
alpha_diversity <- estimate_richness(ms_rare, measures = c("Shannon", "Chao1"))

# Add the 'group' variable to the alpha diversity data
alpha_diversity$group <- sample_data(ms_rare)$group

new_faiths <- data.frame(sample_data(new_faiths))
class(new_faiths)

alpha_diversity$group <- gsub("Control 0", "Healthy", alpha_diversity$group)
alpha_diversity$group <- gsub("Control 1", "Asthma", alpha_diversity$group)
alpha_diversity$group <- gsub("RRMS 0", "MS", alpha_diversity$group)
alpha_diversity$group <- gsub("RRMS 1", "MS + Asthma", alpha_diversity$group)

head(alpha_diversity)

# Create a boxplot for Shannon index (alpha diversity) by 'group'
alpha_diversity_plot <- ggplot(alpha_diversity, aes(x = group, y = Shannon)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (Shannon Index)",
       x = "Group", 
       y = "Shannon Diversity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ # Rotate x-axis labels if needed
  stat_compare_means(method = "kruskal.test", label = "p.signif")
alpha_diversity_plot
ggsave("alpha_diversity_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Shannon ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("alpha_diversity_plot.png", width = 8, height = 6, dpi = 300)


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


