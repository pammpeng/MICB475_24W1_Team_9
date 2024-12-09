library(tidyverse)
library(picante)


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


# Create a data frame with Faith's PD and the group info
faith_pd_df <- data.frame(PD = sample_data(ms_rare)$PD, group = group_info)

# Check the first few rows of the data frame
head(faith_pd_df)

# Perform Kruskal-Wallis test to compare PD across groups
kruskal_test_result <- kruskal.test(PD ~ group, data = faith_pd_df)

# Print the Kruskal-Wallis test result
print(kruskal_test_result)

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ms_rare), aes(group, PD), fill = group) + 
  geom_boxplot() +
  xlab("Group") +
  ylab("Phylogenetic Diversity")+
  theme_minimal() +
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
kruskal_test_result <- kruskal.test(faith_pd ~ group, data = ms_rare))



# Combine two variables (e.g., 'disease_status' and 'treatment') into a new 'group' variable
sample_data(ms_rare)$group <- paste(sample_data(ms_rare)$disease_course, sample_data(ms_rare)$asthma)

sample_variables(ms_rare)

# Calculate alpha diversity
alpha_diversity <- estimate_richness(ms_rare, measures = c("Shannon", "Chao1"))

# Add the 'group' variable to the alpha diversity data
alpha_diversity$group <- sample_data(ms_rare)$group

# Create a boxplot for Shannon index (alpha diversity) by 'group'
alpha_diversity_plot <- ggplot(alpha_diversity, aes(x = group, y = Shannon)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity by Group (Shannon Index)",
       x = "Group", 
       y = "Shannon Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
alpha_diversity_plot

# Run a Kruskal-Wallis test for Shannon diversity
kruskal.test(Shannon ~ group, data = alpha_diversity)

# Save the last plot as a PNG image
ggsave("alpha_diversity_plot.png", width = 8, height = 6, dpi = 300)

