# Code overview

Run every script from the project folder with the `taxonomy` environment active.
Settings for Module 1 are in [`config/pipeline.env`](../config/pipeline.env).

## Module 1 · Read processing

| Script | Role |
| --- | --- |
| `read_processing.py` | **Entry point.** Runs the shell pipeline for every sample, then the analysis. |
| `run_metabarcoding_pipeline.sh` | Per sample: cutadapt + `primer_tags.py` → fastplong → VSEARCH SINTAX. Reuses finished samples. |
| `primer_tags.py` | Reads cutadapt's info file and assigns forward/reverse tags (edit distance, error limit, margin). |
| `sintax_analysis.py` | Joins SINTAX results with tags, contamination screen, per-sample and combined tables. Can be run alone to redo the tables. |
| `sintax_parser.py` | Parses SINTAX output lines into ranks and confidences. |
| `analysis_report.py` | Summary tables and figures at the end of the analysis. |
| `describe_table_headers.py` | Checks every output table and writes `TABLE_HEADERS.md`. |
| `merge_pass_fail.py` | Optional first step: combines MinKNOW `pass`/`fail` subfolders. |
| `reuse_sintax.py` | One-off, for mergedseq5: copies the October 2026 classifications onto the re-cleaned reads, matched by read ID (see [mergedseq5](../data/mergedseq5/README.md)). |
| `sample_folders.py` | Finds `barcodeNN`/`unclassified` folders; stores paths relative to the project. |
| `export_sintax_fasta.py` | Exports reads of a taxon (also used by Module 2). |
| `export_primer_fasta.py` | Exports reads of one forward tag in one barcode folder. |

## Module 2 · Postprocessing

| Script | Role |
| --- | --- |
| `postprocess_sintax.py` | **Entry point.** Export, dominant-taxon selection, clustering, references, protein check, alignment, ClipKIT, IQ-TREE, figure; resumable. |
| `bold_references.py` | Selects named-species records of a taxon from the BOLD FASTA. |
| `genbank_xml.py`, `reference_labels.py` | Read NCBI metadata; label reference IDs with their species. |
| `kcalign_alignment.py` | Runs kc-align with the COI guide and restores the sequence names. |
| `kcalign_cli.py`, `kcalign_parallel.py` | Work around two bugs of the installed kc-align version (thread handling, crash on rejected reads). |
| `plot_phylogeny.py` | Midpoint-rooted tree figure with Baltic. |

## Shared

| Script | Role |
| --- | --- |
| `seqio.py` | FASTA/FASTQ reading and writing. |
| `sequence_headers.py` | Format of local read headers (`folder-forward-reverse-taxon-flag-readid`). |
