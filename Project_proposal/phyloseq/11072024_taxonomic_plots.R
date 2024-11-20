#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(ms_rare, fill="Phylum") 

# Convert to relative abundance
ms_RA <- transform_sample_counts(ms_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
ms_phylum <- tax_glom(ms_RA, taxrank = "Phylum", NArm=FALSE)

plot_bar(ms_phylum, fill="Phylum") + 
  facet_wrap(.~group, scales = "free_x")

gg_taxa <- plot_bar(ms_phylum, fill="Phylum") + 
  facet_wrap(.~group, scales = "free_x")+
  theme_minimal()+
  theme(axis.text.x = element_blank())
gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)


