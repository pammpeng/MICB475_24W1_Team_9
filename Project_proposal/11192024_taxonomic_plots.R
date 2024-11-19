#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(ms_rare, fill="Phylum") 

# Convert to relative abundance
ms_RA <- transform_sample_counts(ms_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
ms_phylum <- tax_glom(ms_RA, taxrank = "Phylum", NArm=FALSE)

# Rename the groups
sample_data(ms_phylum)$group <- gsub("Control_0", "Healthy", sample_data(ms_phylum)$group)
sample_data(ms_phylum)$group <- gsub("Control_1", "Asthma", sample_data(ms_phylum)$group)
sample_data(ms_phylum)$group <- gsub("RRMS_0", "MS", sample_data(ms_phylum)$group)
sample_data(ms_phylum)$group <- gsub("RRMS_1", "MS + Asthma", sample_data(ms_phylum)$group)


plot_bar(ms_phylum, fill="Phylum") + 
  facet_wrap(.~group, scales = "free_x")

gg_taxa <- plot_bar(ms_phylum, fill="Phylum") + 
  facet_wrap(.~group, scales = "free_x")+
  theme_minimal()+
  theme(axis.text.x = element_blank())
gg_taxa

sample_data(ms_rare)$group

head(sample_data(ms_rare))

# Plot again with updated group labels
plot_bar(ms_phylum, fill="Phylum") + 
  facet_wrap(.~group, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)

##Statistical testing for relative abundance

library(vegan)


# Check the Bray-Curtis distance matrix
otu_table_matrix <- as.matrix(otu_table(ms_rare))
bray_dist <- vegdist(t(otu_table_matrix), method = "bray")
bray_dist
# Ensure the dimensions of the Bray-Curtis distance matrix match the number of samples
dim(bray_dist)

# Run PERMANOVA using adonis2
adonis_result <- adonis2(bray_dist ~ group, data = sample_df)

# Print the result
print(adonis_result)
