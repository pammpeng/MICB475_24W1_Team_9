library(tidyverse)
library(picante)
library(ggplot2)
library(dplyr)
library(phyloseq)



install.packages("car")

install.packages("ggpubr", dependencies = TRUE)
library(ggpubr)

#### Alpha diversity ######
plot_richness(Pms_rare)

# -------------------- Preparing Dataframe for alpha diversity plots --------------------

# Extract the grouping variable from sample data
group_info <- sample_data(Pms_rare)$group

#make both disease_course and asthma columns have factors as variables
sample_data(Pms_rare)$disease_course <- as.factor(sample_data(Pms_rare)$disease_course)
sample_data(Pms_rare)$asthma <- as.factor(sample_data(Pms_rare)$asthma)
sample_data(Pms_rare)$group <- as.factor(sample_data(Pms_rare)$group)

# -------------------- Alpha Diversity --------------------
# Calculate alpha diversity
alpha_diversity <- estimate_richness(Pms_rare, measures = c("Shannon", "Chao1", "Observed", "Ace", "Simpson", "Fisher"))

# Add the 'group' variable to the alpha diversity data
alpha_diversity$group <- sample_data(Pms_rare)$group

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

ggsave("shannon_plot.png", width = 8, height = 6, dpi = 300)

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Shannon ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("alpha_diversity_plot.png", width = 8, height = 6, dpi = 300)

