# Sequencing dataset mergedseq5

Raw reads and primer-tracking outputs are organised by sample folder. Cleaned
combined reads are `<folder>/filtered/asinput.fastq.gz`. Classification inputs are
those cleaned reads. Existing SINTAX results live under `data/sintax/<folder>/`.

Read-level analyses use a new `data/sintax_analysis/<run>/` directory and retain
source-folder and PCR-plate provenance. Confidence acceptance and the PCR-plate
prevalence screen both use thresholds of 0.8. No consensus generation or amino-acid
translation/noNUMT filtering is part of this workflow.

See [workflow documentation](../../docs/pipeline_from_fastplong_onwards.md).
