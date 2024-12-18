## Picrust2 QIIME2 Script

# Step 1: Filter features in the feature table based on a minimum frequency
qiime feature-table filter-features \
  --i-table ms_table.qza \  # Input feature table (qiime artifact file), 'ms_table.qza'
  --p-min-frequency 5 \  # Filters out features with a total frequency less than 5 across all samples
  --o-filtered-table feature-frequency-filtered-table.qza  # Output filtered feature table, 'feature-frequency-filtered-table.qza'

# Step 2: Run the PICRUSt2 pipeline to predict functional profiles based on the filtered feature table
qiime picrust2 full-pipeline \
  --i-table feature-frequency-filtered-table.qza \  # Input filtered feature table from Step 1, 'feature-frequency-filtered-table.qza'
  --i-seq ms_rep-seqs.qza \  # Input sequence data (representative sequences), 'ms_rep-seqs.qza'
  --output-dir q2-picrust2_output \  # Output directory where the results will be saved, 'q2-picrust2_output'
  --p-placement-tool sepp \  # Placement method for phylogenetic tree building, using SEPP (SATé-enabled phylogenetic placement)
  --p-hsp-method pic \  # Method for predicting functional genes, using PIC (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
  --p-max-nsti 2 \  # Maximum NSTI (Normalized Strain-Tree Index), the threshold for the reliability of predicted functional genes (higher values indicate less reliable predictions)
  --verbose  # Outputs additional details during the execution of the pipeline for troubleshooting or further insight


#Result said: 2 of 2754 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.

