## Picrust 2 QIIME2 Script


qiime feature-table filter-features \
  --i-table ms_table.qza \
  --p-min-frequency 5 \
  --o-filtered-table feature-frequency-filtered-table.qza

qiime picrust2 full-pipeline \
  --i-table feature-frequency-filtered-table.qza \
  --i-seq ms_rep-seqs.qza \ 
  --output-dir q2-picrust2_output \
  --p-placement-tool sepp \
  --p-hsp-method pic \
  --p-max-nsti 2 \
  --verbose

#2 of 2754 ASVs were above the max NSTI cut-off of 2.0 and were removed from the downstream analyses.

