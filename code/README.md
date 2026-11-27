# Code overview

Run every script from the project folder with the `taxonomy` environment active.
Settings for Module 1 are in [`config/pipeline.env`](../config/pipeline.env).

## Module 1 · Read processing

| Script | Role |
| --- | --- |
| `read_processing.py` | **Entry point.** Runs the shell pipeline for every sample, then the analysis. |
| `run_metabarcoding_pipeline.sh` | Per sample: cutadapt + `primer_tags.py` → fastplong → `cluster_specimens.py` → VSEARCH SINTAX. Reuses finished steps. |
| `primer_tags.py` | Reads cutadapt's info file and assigns forward/reverse tags (edit distance, error limit, margin). |
| `cluster_specimens.py` | Clusters the cleaned reads of each specimen (barcode + forward tag) at 97%: consensus sequences and single reads to classify. |
| `sintax_analysis.py` | One row per cluster with its taxonomy, read count and tags; specimen and species summaries; contamination screen. Can be run alone to redo the tables. |
| `sintax_parser.py` | Reads SINTAX output lines into ranks and confidences. |
| `analysis_report.py` | Species-by-folder table and figures at the end of the analysis. |
| `merge_pass_fail.py` | Optional first step: combines MinKNOW `pass`/`fail` subfolders. |
| `reuse_sintax.py` | One-off, for mergedseq5: reuses the October 2026 classifications of single reads and classifies the other clusters (see [mergedseq5](../data/mergedseq5/README.md)). |
| `sample_folders.py` | Finds `barcodeNN`/`unclassified` folders; stores paths relative to the project; reads settings from `config/pipeline.env`. |
| `export_clusters.py` | Exports the clusters of a taxon and/or forward tag (also used by Module 2); defines the local sequence names. |

## Module 2 · Postprocessing

| Script | Role |
| --- | --- |
| `postprocess_sintax.py` | **Entry point.** Export of clusters, main identification per specimen, references, protein check, alignment, ClipKIT, IQ-TREE, figure; resumable. |
| `bold_references.py` | Selects named-species records of a taxon from the BOLD FASTA. |
| `genbank_xml.py`, `reference_labels.py` | Read NCBI metadata; label reference IDs with their species. |
| `kcalign_alignment.py` | Runs kc-align with the COI guide and restores the sequence names. |
| `kcalign_cli.py`, `kcalign_parallel.py` | Work around two bugs of the installed kc-align version (thread handling, crash on rejected reads). |
| `plot_phylogeny.py` | Midpoint-rooted tree figure with Baltic. |

## Shared

| Script | Role |
| --- | --- |
| `seqio.py` | FASTA/FASTQ/TSV reading and writing; raw read ID of split reads. |
