library(tidyverse)

## reconciling manifest based on filtered metadata
ms_metadata_final<- read_tsv("ms_metadata_final.tsv")
ms_manifest<-read_tsv("ms_manifest.tsv")

# taking only the sample-id column
filtered_metadata_only_samples <- select(ms_metadata_final, `sample-id`)

# joining the filtered sample ids with the manifest file
filtered_manifest <- left_join(filtered_metadata_only_samples, ms_manifest)

# exporting the filtered_manifest to a TSV file
filtered_manifest_filepath <- "filtered_ms_manifest.tsv"
write_tsv(filtered_manifest, filtered_manifest_filepath)
