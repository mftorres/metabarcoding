# Sequencing dataset mergedseq5

Raw reads of sequencing run 5, one folder per ONT barcode (`barcode01`–`barcode24`)
plus `unclassified`.

## Re-analysis with the current pipeline (November 2026)

This dataset was first processed in October 2026 by the previous pipeline version
(primers left on the reads, tags found by mismatch counting). SINTAX had taken weeks,
so it was re-analysed once as follows, instead of classifying again:

1. The old `barcodeNN/filtered/` folders were renamed `filtered_previous_pipeline/`
   (kept), and every sample was cleaned again with the current Module 1
   (`CLEAN_ONLY=1 bash code/run_metabarcoding_pipeline.sh data/mergedseq5`).
2. `python code/reuse_sintax.py mergedseq5 --old-root data/sintax` copied the old
   classifications (`data/sintax/barcodeNN/barcodeNN_sintaxSP.tsv`) onto the new
   cleaned reads, matched by raw read ID, into `data/sintax/mergedseq5/`.
   `reused_classification_report.tsv` there counts, per sample, the reads copied and
   the reads left out: reads the old quality filter had removed, and reads that
   fastplong split into pieces the old run classified differently.
3. `sintax_analysis.py` built the analysis tables from these.

The old classifications were made on reads that still carried primers and tags, so
they can differ slightly from a fresh classification. barcode23 and barcode24 were
never classified. The October 2026 tables (`data/sintax_analysis/mergedseq5_20261021/`)
are in the old format and cannot be used by the current Module 2.
