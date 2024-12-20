library(tidyverse)
library(picante)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(ggsignif)


install.packages("car")

install.packages("ggpubr", dependencies = TRUE)
library(ggpubr)

#### Alpha diversity ######
plot_richness(ms_rare)

# -------------------- Preparing Dataframe for alpha diversity plots --------------------

# Extract the grouping variable from sample data
group_info <- sample_data(ms_rare)$group

#Change group names
sample_data(ms_rare)$group <- gsub("Control_0", "Healthy", sample_data(ms_rare)$group)
sample_data(ms_rare)$group <- gsub("Control_1", "Asthma", sample_data(ms_rare)$group)
sample_data(ms_rare)$group <- gsub("RRMS_0", "MS", sample_data(ms_rare)$group)
sample_data(ms_rare)$group <- gsub("RRMS_1", "MS + Asthma", sample_data(ms_rare)$group)
#view updated names
sample_data(ms_rare)

#make both disease_course and asthma columns have factors as variables
sample_data(ms_rare)$disease_course <- as.factor(sample_data(ms_rare)$disease_course)
sample_data(ms_rare)$asthma <- as.factor(sample_data(ms_rare)$asthma)
sample_data(ms_rare)$group <- as.factor(sample_data(ms_rare)$group)

# -------------------- Faiths Phylogenetic Diversity --------------------

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(ms_rare)), phy_tree(ms_rare),
                 include.root=F) 
?pd

# add PD to metadata table
sample_data(ms_rare)$PD <- phylo_dist$PD
sample_data(ms_rare)
sample_data(ms_rare)$group <- paste(sample_data(ms_rare)$disease_course, sample_data(ms_rare)$asthma, sep = "_")


# Create a data frame with Faith's PD and the group info
faith_pd_df <- data.frame(PD = sample_data(ms_rare)$PD, group = group_info)

# Check the first few rows of the data frame
head(faith_pd_df)

# Print the Kruskal-Wallis test result
print(kruskal_test_result)

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ms_rare), aes(group, PD), fill = group) + 
  geom_boxplot() +
  xlab("Group") +
  ylab("Phylogenetic Diversity")+
  theme_classic() 
# view plot
plot.pd

# Save the last plot as a PNG image
ggsave("plot_pd.png", width = 8, height = 6, dpi = 300)


# -------------------- Alpha Diversity --------------------
# Calculate alpha diversity
alpha_diversity <- estimate_richness(ms_rare, measures = c("Shannon", "Chao1", "Observed", "Ace", "Simpson", "Fisher"))

# Add the 'group' variable to the alpha diversity data
alpha_diversity$group <- sample_data(ms_rare)$group

#rename groups
alpha_diversity$group <- gsub("Control_0", "Healthy", alpha_diversity$group)
alpha_diversity$group <- gsub("Control_1", "Asthma", alpha_diversity$group)
alpha_diversity$group <- gsub("RRMS_0", "MS", alpha_diversity$group)
alpha_diversity$group <- gsub("RRMS_1", "MS + Asthma", alpha_diversity$group)
#check group names
head(alpha_diversity)

alpha_diversity

## Shannon Index Boxplot ##

# Create a boxplot for Shannon index (alpha diversity) by 'group'
shannon_plot <- ggplot(alpha_diversity, aes(x = group, y = Shannon)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))+
  theme(axis.title.x = element_blank())+
  theme(plot.title = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 14))
shannon_plot

shannon_plot <- ggplot(alpha_diversity, aes(x = group, y = Shannon, color = group)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
  theme_classic()+
  theme(axis.title.x = element_blank())+
  labs(y = "Shannon Diversity Index")+
  theme(axis.title.y = element_text(size = 16))+
  theme(axis.text.x = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size = 16))+
  ylim(1, 5)
  geom_signif(
    comparisons = list(c("group", "shannon"), c("MS", "Healthy"), c("Asthma", "MS + Asthma")),
    map_signif_level = TRUE,
    test = "kruskal.test",
    textsize = 6,
    y_position = 3)
    
shannon_plot

# Install and load the dunn.test package for post-hoc pairwise tests
install.packages("dunn.test")
library(dunn.test)

# Perform pairwise comparisons if Kruskal-Wallis is significant
dunn_result <- dunn.test(alpha_diversity$Shannon, alpha_diversity$group, 
                         kw = TRUE, label = TRUE)

# View the Dunn test result
print(dunn_result)



shannon_plot

ggsave("shannon_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Shannon ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("alpha_diversity_plot.png", width = 8, height = 6, dpi = 300)

## Chao1 Index Boxplot ##
chao1_plot <- ggplot(alpha_diversity, aes(x = group, y = Chao1)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (Chao1 Index)",
       y = "Alpha Diversity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x = element_blank())
chao1_plot

ggsave("chao1_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Chao1 ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("chao1_plot.png", width = 8, height = 6, dpi = 300)

## Observed Index Boxplot ##
observed_plot <- ggplot(alpha_diversity, aes(x = group, y = Observed)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (Observed Index)",
       y = "Alpha Diversity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x = element_blank())
observed_plot

ggsave("observed_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Observed ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("observed_plot.png", width = 8, height = 6, dpi = 300)

## Ace Index Boxplot ##
ace_plot <- ggplot(alpha_diversity, aes(x = group, y = ACE)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (ACE Index)",
       y = "Alpha Diversity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x = element_blank())
ace_plot

ggsave("ace_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(ACE ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("ace_plot.png", width = 8, height = 6, dpi = 300)

## Fisher Index Boxplot ##
fisher_plot <- ggplot(alpha_diversity, aes(x = group, y = Fisher)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (Fisher Index)",
       y = "Alpha Diversity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x = element_blank())
fisher_plot

ggsave("fisher_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Fisher ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("fisher_plot.png", width = 8, height = 6, dpi = 300)

## Simpson Index Boxplot ##
simpson_plot <- ggplot(alpha_diversity, aes(x = group, y = Simpson)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (Simpson Index)",
       y = "Alpha Diversity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
simpson_plot

ggsave("simpson_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Simpson ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("simpson_plot.png", width = 8, height = 6, dpi = 300)

