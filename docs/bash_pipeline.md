# Bash Metabarcoding Pipeline

This is a Snakemake-free fallback pipeline. It runs the same stages:

1. `fastplong` on each `*.fastq.gz` file in a barcode folder
2. concatenate filtered FASTQ files into `filtered/combined.fastq.gz`
3. `code/amplicon_sorter.py`
4. `kraken2` on amplicon-sorter FASTA outputs

Edit defaults in:

```bash
config/pipeline.env
```

Run one folder:

```bash
conda activate saproxdna
bash code/run_metabarcoding_pipeline.sh data/merged/barcode58
```

Run three folders with two folders active at a time:

```bash
conda activate saproxdna
MAX_JOBS=2 bash code/run_metabarcoding_pipeline.sh \
  data/merged/barcode58 \
  data/merged/barcode60 \
  data/merged/barcode61
```

Run three folders at once:

```bash
MAX_JOBS=3 bash code/run_metabarcoding_pipeline.sh \
  data/merged/barcode58 \
  data/merged/barcode60 \
  data/merged/barcode61
```

Override parameters without editing the config file:

```bash
MAX_JOBS=3 \
AMPLICON_THREADS=12 \
KRAKEN_THREADS=10 \
AMPLICON_SG=99.0 \
AMPLICON_SS=99.0 \
AMPLICON_SC=99.0 \
KRAKEN_DB=/mnt/d/kraken2DB/k2_core_nt_2025 \
bash code/run_metabarcoding_pipeline.sh data/merged/barcode58 data/merged/barcode60
```

The script writes `.done` files in each folder's `filtered/` directory and skips completed steps. To force rerun:

```bash
FORCE=1 bash code/run_metabarcoding_pipeline.sh data/merged/barcode58
```

To rerun Kraken2 only, without rerunning filtering, combining, or amplicon-sorter:

```bash
FORCE_KRAKEN=1 bash code/run_metabarcoding_pipeline.sh data/merged/barcode58
```

Kraken2 target FASTA files are concatenated into one file per sample folder before classification:

```bash
data/merged/barcodeXX/filtered/combined/kraken2_combined_input.fasta
```

The matching Kraken2 output and report are written as `kraken2_combined_input.fasta.kraken` and `kraken2_combined_input.fasta.report`.

If amplicon-sorter outputs already exist but `amplicon_sorter.done` is missing, the script now detects the output files and recreates the `.done` marker instead of rerunning amplicon-sorter.

Logs are written under:

```bash
data/merged/barcodeXX/filtered/logs/
```
