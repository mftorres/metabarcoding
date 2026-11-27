# Sequencing dataset mergedseq5

Raw reads of sequencing run 5, one folder per ONT barcode (`barcode01`–`barcode24`)
plus `unclassified`.

## Re-analysis with the current pipeline (November 2026)

This dataset was first processed in October 2026 by the previous pipeline version
(primers left on the reads, tags found by mismatch counting, every read classified).
SINTAX had taken weeks, so it was re-analysed once as follows instead of classifying
everything again:

1. The old `barcodeNN/filtered/` folders were renamed `filtered_previous_pipeline/`
   (kept). Every sample was cleaned and clustered again with the current Module 1,
   without classification:
   `CLEAN_ONLY=1 bash code/run_metabarcoding_pipeline.sh data/mergedseq5`.
2. `python code/reuse_sintax.py mergedseq5 --old-root data/sintax --db data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.udb`
   wrote one classification per cluster to `data/sintax/mergedseq5/`. A cluster of
   one read is the same read that was classified in October
   (`data/sintax/barcodeNN/barcodeNN_sintaxSP.tsv`), so its October classification
   was copied, matched by raw read ID. Consensus sequences, and single reads without
   a usable October result, were classified with the same full MIDORI database.
   `reused_classification_report.tsv` counts both per sample.
3. barcode23 and barcode24 had never been classified; they were classified normally
   with the full MIDORI database (`SINTAX_SUFFIX=sintaxSP`).
4. `python code/sintax_analysis.py mergedseq5 --suffix sintaxSP --output-dir data/read_processing/mergedseq5_clusters/analysis`
   built the analysis tables. Use this folder as `--analysis-dir` in Module 2.

The October classifications were made on reads that still carried primers and tags,
so they can differ slightly from a fresh classification.

Older results, kept but not usable by the current Module 2:

| Folder | What it is |
| --- | --- |
| `data/sintax_analysis/mergedseq5_20261021/` | October 2026 tables, old format. |
| `data/sintax/mergedseq5_readlevel_20261123/` | October classifications copied onto the re-cleaned reads, one per read (23 November 2026). |
| `data/read_processing/mergedseq5_reused/` | Analysis built from those, one row per read, before clustering moved to Module 1. |
