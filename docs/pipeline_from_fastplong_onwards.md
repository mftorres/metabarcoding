# Pipeline from fastplong onwards

This document records the current workflow from raw per-barcode FASTQ files through fastplong filtering, amplicon_sorter grouping, Kraken2 classification, primer/UMI tracking, read-level tracking-table construction, and final evaluation-table export.

The examples below assume the project root is:

```bash
cd /mnt/c/Users/mftor/Documents/metabarcoding
```

## 0. Required starting layout

Raw reads should be directly inside one folder per sequencing barcode:

```text
data/
  mergedseq5/
    barcode01/
      *.fastq.gz
    barcode02/
      *.fastq.gz
    ...
```

If files arrive under `pass/barcodeXX` or `fail/barcodeXX`, move or copy the relevant FASTQ files into `data/mergedseqX/barcodeXX/` before starting this pipeline.

## 1. Configure the shell pipeline

Default shell-pipeline settings are in:

```text
config/pipeline.env
```

Important variables:

| Variable | Purpose |
|---|---|
| `ADAPTER_FASTA` | ONT barcode/adapter FASTA used by fastplong. |
| `FASTPLONG_UNQUALIFIED_PERCENT_LIMIT` | fastplong maximum unqualified-base percentage. |
| `FASTPLONG_QUALIFIED_QUALITY_PHRED` | fastplong qualified-base quality threshold. |
| `FASTPLONG_MEAN_QUAL` | fastplong mean-quality threshold. |
| `FASTPLONG_LENGTH_REQUIRED` | fastplong minimum read length. Current default is `300`. |
| `FASTPLONG_LENGTH_LIMIT` | fastplong maximum read length. Current default is `1500`. |
| `AMPLICON_SORTER_SCRIPT` | Local amplicon_sorter script. This should be `code/amplicon_sorter_mod.py` because it writes the read map needed downstream. |
| `AMPLICON_READ_MAP_SCRIPT` | Script that builds the Kraken-input-to-ONT-read map. |
| `AMPLICON_THREADS` | Threads used by amplicon_sorter. |
| `AMPLICON_MINLENGTH` | Minimum length passed to amplicon_sorter. |
| `AMPLICON_MAXLENGTH` | Maximum length passed to amplicon_sorter. |
| `AMPLICON_SG` | amplicon_sorter similar-genes threshold. |
| `AMPLICON_SS` | amplicon_sorter similar-species-groups threshold. |
| `AMPLICON_SC` | amplicon_sorter similar-consensus threshold. |
| `AMPLICON_TMP_ROOT` | Temporary directory for amplicon_sorter runs. |
| `KRAKEN_THREADS` | Threads used by Kraken2. |
| `KRAKEN_DB` | Kraken2 database directory. |
| `KRAKEN_CONFIDENCE` | Kraken2 confidence threshold. |
| `KRAKEN_MEMORY_MAPPING` | Whether to use Kraken2 memory mapping. |
| `KRAKEN_USE_DOCKER` | Set to `1` to run Kraken2 through Docker instead of local `kraken2`. |
| `KRAKEN_DOCKER_IMAGE` | Docker image used when `KRAKEN_USE_DOCKER=1`. |
| `MAX_JOBS` | Number of barcode folders processed concurrently by shell wrappers. |
| `FORCE` | Set to `1` to rerun completed steps. |
| `FORCE_KRAKEN` | Set to `1` to rerun only Kraken-related steps while preserving earlier outputs. |
| `PRIMER_UMI_FASTA` | Primer/UMI FASTA for seqkit primer tracking. |
| `PRIMER_UMI_MAX_MISMATCH` | Maximum mismatches allowed by `seqkit locate`. Current default is `10`; substantially larger values can generate extremely large weak-match tables. |
| `PRIMER_UMI_SEQKIT_THREADS` | Threads used by seqkit primer tracking. |
| `NONUMT_MIN_AA_LENGTH` | Minimum amino-acid length kept by the consensus noNUMT filter. |
| `NONUMT_TRANSLATION_TABLE` | Genetic code table used by `seqkit translate` for the consensus noNUMT filter. |
| `NONUMT_SEQKIT_THREADS` | Threads used by `seqkit translate` for the consensus noNUMT filter. |
| `SEQKIT_CONDA_ENV` | Conda environment used by noNUMT filtering when `seqkit` is not available on `PATH`. Current default is `umitracker`. |
| `SEQKIT_ENV` | Conda environment used by `code/run_full_downstream_pipeline.sh` for primer/UMI tracking. Current default is `umitracker`. |
| `SHELL_ENV` | Conda environment used by `code/run_full_downstream_pipeline.sh` for shell/Python helper scripts that do not require seqkit. Current default is `saproxdna`. |
| `NOTEBOOK_ENV` | Conda environment used by `code/run_full_downstream_pipeline.sh` for notebooks and pandas-heavy Python scripts. Current default is `saproxdna-evaluate`. |

Activate the environment that has fastplong, amplicon_sorter dependencies, and Kraken2 for the shell pipeline. `seqkit` is currently in the `umitracker` conda environment; `code/run_metabarcoding_pipeline.sh` will use `seqkit` from the active environment if available, otherwise it falls back to `conda run -n "$SEQKIT_CONDA_ENV" seqkit` for the noNUMT step. In the current work, notebook evaluation uses the `saproxdna-evaluate` environment.

The full downstream wrapper switches conda environments by step. `code/run_full_downstream_pipeline.sh` uses `SHELL_ENV=saproxdna` for the fastplong/amplicon_sorter/noNUMT/Kraken2 pipeline and shell helpers, `SEQKIT_ENV=umitracker` for primer/UMI tracking, and `NOTEBOOK_ENV=saproxdna-evaluate` for notebooks by default. It first runs `code/run_metabarcoding_pipeline.sh` for the selected barcode range, so deleted or missing `filtered/` outputs are rebuilt while valid completed outputs are skipped. Override these variables only if packages are moved between environments.

## 2. Run fastplong, amplicon_sorter, Kraken2, and Kraken read-map generation

Script:

```text
code/run_metabarcoding_pipeline.sh
```

This script runs, per barcode:

1. fastplong quality/filter/adapter trimming, including `--length_required "$FASTPLONG_LENGTH_REQUIRED"` and `--length_limit "$FASTPLONG_LENGTH_LIMIT"`.
2. Combine filtered FASTQ files into `filtered/asinput.fastq.gz`.
3. amplicon_sorter using `code/amplicon_sorter_mod.py`.
4. Translate consensus sequences with `seqkit` and write cleaned nucleotide noNUMT consensus FASTA files.
5. Build a read map for the single cleaned `*consensussequences_noNUMT_*filter.fasta` file.
6. Kraken2 classification of the cleaned noNUMT consensus FASTA.

If amplicon_sorter finishes successfully but assigns no reads to a final
consensus group, the barcode is recorded as a valid no-consensus result:

```text
data/mergedseq5/barcodeXX/filtered/no_consensus.done
```

The pipeline retains the amplicon_sorter results/read map, skips noNUMT output
generation and Kraken2 classification for that barcode, and continues with the
remaining barcodes. Read/primer tracking is still generated with blank
consensus taxIDs. The per-barcode taxonomy evaluation notebook is skipped
because there is no consensus taxonomy to evaluate, and the phylogenetic
exporter writes an empty per-barcode FASTA.

Run one barcode:

```bash
conda activate saproxdna
bash code/run_metabarcoding_pipeline.sh data/mergedseq5/barcode01
```

Run all barcodes in a folder:

```bash
conda activate saproxdna
MAX_JOBS=2 bash code/run_metabarcoding_pipeline.sh data/mergedseq5/barcode{01..24}
```

Rerun everything for one barcode:

```bash
FORCE=1 bash code/run_metabarcoding_pipeline.sh data/mergedseq5/barcode01
```

Rerun Kraken2/read-map parts while keeping filtering and amplicon_sorter outputs:

```bash
FORCE_KRAKEN=1 bash code/run_metabarcoding_pipeline.sh data/mergedseq5/barcode01
```

Main per-barcode outputs:

```text
data/mergedseq5/barcodeXX/filtered/
  *_filtered.fastq.gz
  *_report.html
  asinput.fastq.gz
  filtering.done
  combined.done
  amplicon_sorter.done
  kraken2.done
  consensusfile.fasta
  results.csv
  logs/
    *.log
    amplicon_sorter.stdout.log
    amplicon_sorter.stderr.log
    amplicon_sorter.log
    kraken2.log
  combined/
    amplicon_sorter_read_map.tsv
    results.txt
    asinput_*_*.fasta
    asinput_*_unique.fasta
    asinput_consensussequences.fasta
    asinput_consensussequences_aa_noNUMT_200filter.fasta
    asinput_consensussequences_noNUMT_200filter.fasta
    asinput_consensussequences_noNUMT_200filter.read_map.tsv
    asinput_consensussequences_noNUMT_200filter.fasta.kraken
    asinput_consensussequences_noNUMT_200filter.fasta.report
```

### Important read-map outputs

`filtered/combined/amplicon_sorter_read_map.tsv` links amplicon_sorter internal read numbers to ONT read IDs.

Before noNUMT filtering and Kraken2, the pipeline prefixes every consensus FASTA
header with its barcode (for example, `barcode01_consensus_asinput_0_0(1036)`).
This is idempotent, so rerunning the pipeline does not add the prefix twice. It
prevents identical amplicon_sorter consensus IDs from different barcodes being
collapsed or assigned the wrong Kraken2 result when files are combined. The
barcode-prefixed ID is retained in Kraken outputs and read maps; parsers accept
both the new prefixed form and older unprefixed records.

The current `run_metabarcoding_pipeline.sh` invocation runs Kraken2 separately
for each barcode. If consensus FASTAs are later concatenated for a merged Kraken2
run, retain these prefixes; they are the uniqueness boundary between barcodes.

| Column | Meaning |
|---|---|
| `amplicon_sorter_read_number` | Numeric read identifier assigned internally by amplicon_sorter. |
| `ont_read_id` | Original or processed ONT read ID. Split reads can be named `split-by-adapter-left-<raw_id>` or `split-by-adapter-right-<raw_id>`. |
| `read_length` | Length of the sequence passed into amplicon_sorter. |

`filtered/combined/asinput_consensussequences_noNUMT_200filter.read_map.tsv` links every Kraken input record to amplicon_sorter context.

| Column | Meaning |
|---|---|
| `merge_folder` | Parent sequencing merge/run folder, for example `mergedseq5`. |
| `barcode` | Barcode folder, for example `barcode01`. |
| `kraken_read_id` | Consensus FASTA record ID sent to Kraken2. |
| `amplicon_sorter_group` | Group associated with the consensus record, for example `1_0`. |
| `kraken_sequence_length` | Length of the exact sequence sent to Kraken2. |

## 3. Check amplicon_sorter FASTA sequence lengths

Script:

```text
code/length_check_amplicon_sorter.sh
```

Run this after amplicon_sorter and Kraken2 have completed for the barcode folders. It can be run after the full `code/run_metabarcoding_pipeline.sh` step because that pipeline includes amplicon_sorter, Kraken2, and the Kraken read map.

The script scans all barcodes in a merged-sequence folder:

```text
data/mergedseq5/barcode*/filtered/combined/asinput_*_*.fasta
```

It also screens the per-barcode Kraken2 mapping/classification files:

```text
data/mergedseq5/barcode*/filtered/combined/asinput_consensussequences_noNUMT_200filter.read_map.tsv
data/mergedseq5/barcode*/filtered/combined/asinput_consensussequences_noNUMT_200filter.fasta.kraken
data/mergedseq5/barcode*/filtered/combined/asinput_consensussequences_noNUMT_200filter.fasta.report
```

Kraken2 taxid assignment is attached as follows:

- For grouped files such as `asinput_1_0.fasta`, the script uses the corresponding consensus sequence in `asinput_consensussequences_noNUMT_200filter.fasta.kraken`. The link is made through `asinput_consensussequences_noNUMT_200filter.read_map.tsv`, where consensus read IDs follow the `consensus_combinedrawreads_<cluster>_<group>(...)` pattern.
- Unique-read FASTA files are not sent to Kraken2 in the current consensus-only pipeline, so they do not receive direct Kraken taxids in this table.
- Taxon names are taken from `asinput_consensussequences_noNUMT_200filter.fasta.report`.

Run with defaults:

```bash
code/length_check_amplicon_sorter.sh
```

The default input folder is `data/mergedseq5`, and the default output table is:

```text
data/mergedseq5/mergedseq5_amplicon_sorter_length_check.tsv
```

Run with explicit input and output paths:

```bash
code/length_check_amplicon_sorter.sh \
  data/mergedseq5 \
  data/mergedseq5/mergedseq5_amplicon_sorter_length_check.tsv
```

The output table is tab-delimited and contains one row per `asinput_*_*.fasta` file.

| Column | Meaning |
|---|---|
| `barcode` | Barcode folder containing the FASTA file. |
| `asinput_file` | FASTA filename, for example `asinput_1_0.fasta` or `asinput_0_unique.fasta`. |
| `cluster_group` | Parsed cluster/group label without the `asinput_` prefix or `.fasta` suffix, for example `1_0` or `0_unique`. |
| `cluster_number` | First parsed value after `asinput_`. |
| `group_number` | Second parsed value after `asinput_`. For unique files, this is `unique`. |
| `sequence_count` | Number of FASTA records in the file. |
| `average_length` | Mean sequence length in bases. |
| `length_standard_deviation` | Sample standard deviation of sequence lengths in bases. Files with one sequence have standard deviation `0`. |
| `kraken_record_count` | Number of consensus Kraken records associated with this FASTA file. Grouped files usually have one consensus record; unique files have `0` because they are not sent to Kraken2. |
| `kraken_classified_count` | Number of associated Kraken records classified by Kraken2. |
| `kraken_unclassified_count` | Number of associated Kraken records unclassified by Kraken2. |
| `kraken_taxids` | Semicolon-separated Kraken taxids associated with the file. Taxid `0` means unclassified. |
| `kraken_taxon_names` | Semicolon-separated taxon names corresponding to `kraken_taxids`, in the same order. Names are taken from `asinput_consensussequences_noNUMT_200filter.fasta.report`. |
| `kraken_taxid_counts` | Semicolon-separated `taxid:count` pairs for the associated Kraken records. |
| `file_path` | Full path to the FASTA file used for the calculation. |

## 4. noNUMT consensus filtering before Kraken2

This step is part of the main shell pipeline:

```text
code/run_metabarcoding_pipeline.sh
```

It runs after amplicon_sorter has produced per-barcode consensus FASTA files:

```text
data/mergedseq5/barcode*/filtered/combined/asinput_consensussequences.fasta
```

The pipeline uses each amplicon_sorter consensus FASTA as input. It translates each consensus sequence in all six frames, appends the frame to the translated amino-acid header, trims terminal stop codons, and keeps translated amino-acid sequences at or above `NONUMT_MIN_AA_LENGTH`.

The translated amino-acid FASTA is used only as a diagnostic/filtering intermediate. The pipeline then strips any `_frame=#` or `_frame=-#` text from passing translated headers, deduplicates the passing consensus IDs, and extracts those original nucleotide consensus sequences into a cleaned noNUMT nucleotide FASTA. Kraken2 uses this cleaned nucleotide FASTA, not the amino-acid FASTA.

Useful variables from `config/pipeline.env`:

| Variable | Default | Meaning |
|---|---:|---|
| `NONUMT_MIN_AA_LENGTH` | `200` | Minimum amino-acid length retained by `seqkit seq --min-len`. |
| `NONUMT_TRANSLATION_TABLE` | `5` | Genetic code table passed to `seqkit translate -T`. |
| `NONUMT_SEQKIT_THREADS` | `6` | Threads passed to `seqkit translate -j`. |
| `SEQKIT_CONDA_ENV` | `umitracker` | Conda environment used when `seqkit` is not available in the active shell. |
| `FORCE` | `0` | Set to `1` to overwrite existing noNUMT outputs. |

For each input:

```text
asinput_consensussequences.fasta
```

the pipeline writes two files in the same folder:

```text
asinput_consensussequences_aa_noNUMT_200filter.fasta
asinput_consensussequences_noNUMT_200filter.fasta
```

The `aa` file keeps translated amino-acid sequences and may contain `_frame=...` in the headers. The non-`aa` file contains nucleotide consensus sequences with cleaned headers and is the file used as Kraken2 input.

## 5. Run primer/UMI tracking with seqkit

Script:

```text
code/run_primer_umi_tracking.sh
```

This script runs `seqkit locate` on raw FASTQ files and writes one primer-hit table per barcode.

Run one barcode:

```bash
conda activate saproxdna
bash code/run_primer_umi_tracking.sh data/mergedseq5/barcode01
```

Run all barcode folders under one parent:

```bash
conda activate saproxdna
MAX_JOBS=4 bash code/run_primer_umi_tracking.sh data/mergedseq5
```

Per-barcode outputs:

```text
data/mergedseq5/barcodeXX/
  rawreads_barcodeXX.txt
  barcodeXX_primer_hits.tsv
  barcodeXX_primer_hits.log
```

`barcodeXX_primer_hits.tsv` columns:

| Column | Meaning |
|---|---|
| `seqID` | Raw read identifier searched by seqkit. |
| `patternName` | Primer/UMI pattern name detected. |
| `pattern` | Primer/UMI sequence pattern. |
| `strand` | Strand of the detected match. |
| `start` | Start coordinate of the match. |
| `end` | End coordinate of the match. |
| `matched` | Sequence segment matching the primer/UMI pattern. |

## 6. Parse primer/UMI hits into best-hit tables

Notebook:

```text
code/primer_umi_seqkit_parser.ipynb
```

The notebook itself reads each `barcodeXX_primer_hits.tsv` in chunks, keeps only
the best forward and reverse hit for each read in memory, and writes the outputs
incrementally. This prevents the notebook kernel from attempting to load tens
of gigabytes of raw SeqKit hits at once. The all-hit table is an on-disk audit
output and is not loaded for visualisation. Set `WRITE_ALL_HITS = False` in the
notebook when that very large audit table is not needed; adjust `CHUNK_SIZE`
there if memory use needs to be tuned.

This notebook parses per-barcode `barcodeXX_primer_hits.tsv` files and writes combined primer/UMI tracking outputs under:

```text
data/primertrack/
  mergedseqX_seqkit_primer_umi_all_hits_with_metadata_YYYYMMDD.txt
  mergedseqX_seqkit_primer_umi_best_hits_YYYYMMDD.txt
  mergedseqX_seqkit_primer_umi_best_hit_counts_by_barcode_YYYYMMDD.txt
```

The downstream read-tracking script needs the best-hit file:

```text
data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_20260630.txt
```

The best-hit file is tab-delimited. Important columns used downstream are:

| Column | Meaning |
|---|---|
| `merged_folder` | Merge/run folder, for example `mergedseq5`. |
| `barcode_folder` | Barcode folder, for example `barcode01`. |
| `seqID` | Raw ONT read ID. |
| `forward_patternName` | Best forward primer/UMI hit name. |
| `forward_strand` | Strand of the best forward hit. |
| `forward_start` | Start coordinate of the best forward hit. |
| `forward_end` | End coordinate of the best forward hit. |
| `forward_mismatch_count` | Mismatch count for the best forward hit. |
| `forward_primer_identity` | Identity score for the best forward hit. |
| `reverse_patternName` | Best reverse primer/UMI hit name. |
| `reverse_strand` | Strand of the best reverse hit. |
| `reverse_start` | Start coordinate of the best reverse hit. |
| `reverse_end` | End coordinate of the best reverse hit. |
| `reverse_mismatch_count` | Mismatch count for the best reverse hit. |
| `reverse_primer_identity` | Identity score for the best reverse hit. |

## 7. Check amplicon_sorter primer diversity

Script:

```text
code/primer_diversity_amplicon_sorter.sh
```

Run this after amplicon_sorter has produced the per-barcode `filtered/combined/asinput_*_*.fasta` files and after `code/primer_umi_seqkit_parser.ipynb` has produced the combined seqkit best-hit table.

The script uses:

```text
data/mergedseq5/barcode*/filtered/combined/asinput_*_*.fasta
data/mergedseq5/barcode*/filtered/combined/amplicon_sorter_read_map.tsv
data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_YYYYMMDD.txt
```

The script translates amplicon_sorter FASTA record IDs back to original ONT read IDs using `amplicon_sorter_read_map.tsv`, then joins those ONT read IDs to the seqkit best-hit table. This gives the detected forward- and reverse-primer diversity within each grouped or unique amplicon_sorter FASTA file.

Run with defaults:

```bash
code/primer_diversity_amplicon_sorter.sh
```

The default command automatically uses the latest matching best-hit table under `data/primertrack/`. To force a specific dated best-hit table:

```bash
code/primer_diversity_amplicon_sorter.sh \
  data/mergedseq5 \
  data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_YYYYMMDD.txt
```

The default output table is:

```text
data/mergedseq5/mergedseq5_amplicon_sorter_primer_diversity.tsv
```

The output table is tab-delimited and contains one row per `asinput_*_*.fasta` file.

| Column | Meaning |
|---|---|
| `barcode` | Barcode folder containing the FASTA file. |
| `asinput_file` | FASTA filename, for example `asinput_1_0.fasta` or `asinput_0_unique.fasta`. |
| `cluster_group` | Parsed cluster/group label without the `asinput_` prefix or `.fasta` suffix. |
| `cluster_number` | First parsed value after `asinput_`. |
| `group_number` | Second parsed value after `asinput_`. For unique files, this is `unique`. |
| `sequence_count` | Number of FASTA records in the asinput file. |
| `raw_read_count` | Number of FASTA records successfully traced back to original ONT read IDs. |
| `forward_detected_read_count` | Number of traced raw reads with a detected forward-primer best hit in the seqkit best-hit table. |
| `reverse_detected_read_count` | Number of traced raw reads with a detected reverse-primer best hit in the seqkit best-hit table. |
| `unique_forward_primer_count` | Number of distinct detected forward-primer names in the asinput file. This is the forward-primer diversity measure. |
| `unique_reverse_primer_count` | Number of distinct detected reverse-primer names in the asinput file. This is the reverse-primer diversity measure. |
| `forward_primers_detected` | Semicolon-separated list of distinct detected forward-primer names. |
| `reverse_primers_detected` | Semicolon-separated list of distinct detected reverse-primer names. |
| `forward_primer_counts` | Semicolon-separated `forward_primer:read_count` pairs. |
| `reverse_primer_counts` | Semicolon-separated `reverse_primer:read_count` pairs. |
| `file_path` | Full path to the FASTA file used for the calculation. |

## 8. Parse Kraken2 outputs

Notebook:

```text
code/kraken_parser.ipynb
```

This notebook reads Kraken2 outputs and noNUMT consensus read-map files, attaches barcode/run metadata, and exports parsed Kraken/taxonomy tables under:

```text
data/kraken_parsed/
  mergedseqX_all_kraken_reads_YYYYMMDD.csv
  mergedseqX_<target>_kraken_reads_YYYYMMDD.csv
  mergedseqX_<target>_taxids_YYYYMMDD.csv
  mergedseqX_<target>_read_counts_by_taxon_YYYYMMDD.csv
  mergedseqX_<target>_read_counts_by_barcode_YYYYMMDD.csv
  mergedseqX_kraken_files_missing_reports_YYYYMMDD.csv
```

These files currently retain a `.csv` extension but are tab-delimited. Load them with:

```python
pd.read_csv(path, sep='\t')
```

Important: these Kraken parser outputs are dated. Every time `code/kraken_parser.ipynb` is rerun, it creates new filenames using the current date, for example `mergedseq5_all_kraken_reads_YYYYMMDD.csv` and `mergedseq5_insecta_kraken_reads_YYYYMMDD.csv`. The later steps do not automatically know which dated parser output is the intended one. After creating a new parser output, update every later command or code cell that points to these files, especially the `--insecta-kraken-reads` argument passed to `code/build_barcode_read_tracking_table.py`.

The read-tracking script currently uses the Insecta reads file:

```text
data/kraken_parsed/mergedseq5_insecta_kraken_reads_20260630.csv
```

Required columns in that file:

| Column | Meaning |
|---|---|
| `readID` | Kraken input consensus record ID. |
| `merge_folder` | Merge/run folder. |
| `barcode` | Barcode folder. |

Additional taxonomy columns are produced by the parser and used for reporting, but `build_barcode_read_tracking_table.py` only needs the three columns above to intersect reads/groups with a target clade such as Insecta.

## 9. Build per-read primer/group/taxid tracking tables

Script:

```text
code/build_barcode_read_tracking_table.py
```

This script combines:

- raw ONT read IDs and raw read lengths;
- seqkit best forward/reverse primer hits;
- amplicon_sorter processed read IDs and split-read fragments;
- amplicon_sorter subgroup membership;
- Kraken2 group-consensus taxids propagated to reads assigned to each amplicon_sorter group;
- the Insecta subset from `kraken_parser.ipynb`.

Run one barcode:

```bash
python code/build_barcode_read_tracking_table.py \
  --barcode-dir data/mergedseq5/barcode01 \
  --primer-best-hits data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_20260630.txt \
  --insecta-kraken-reads data/kraken_parsed/mergedseq5_insecta_kraken_reads_20260630.csv \
  --output data/mergedseq5/barcode01/barcode01_read_primer_group_taxid.tsv
```

If you regenerated the parsed Kraken tables, replace the dated file in `--insecta-kraken-reads` with the newest intended `mergedseq5_insecta_kraken_reads_YYYYMMDD.csv`. Using an old dated file will silently produce stale or incomplete Insecta-tracking tables for any barcodes missing from that old parser output.

Run all barcodes in a folder:

```bash
for barcode_dir in data/mergedseq5/barcode{01..24}; do
  barcode="$(basename "$barcode_dir")"
  echo "Building read tracking for ${barcode}"

  python code/build_barcode_read_tracking_table.py \
    --barcode-dir "$barcode_dir" \
    --primer-best-hits data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_20260630.txt \
    --insecta-kraken-reads data/kraken_parsed/mergedseq5_insecta_kraken_reads_20260630.csv \
    --output "$barcode_dir/${barcode}_read_primer_group_taxid.tsv"
done
```

The script writes four outputs per barcode:

```text
barcodeXX_read_primer_group_taxid.tsv.gz
barcodeXX_read_primer_group_taxid_reduced.tsv
barcodeXX_insecta_read_primer_group_taxid.tsv.gz
barcodeXX_insecta_read_primer_group_taxid_reduced.tsv
```

The full tables are compressed after the reduced tables are written. The reduced tables remain uncompressed because the evaluation notebook loads them faster.

Full read-tracking table columns:

| Column | Meaning |
|---|---|
| `merge_folder` | Sequencing merge/run folder, for example `mergedseq5`. |
| `barcode` | Barcode folder, for example `barcode01`. |
| `raw_ont_read_id` | Original ONT read UUID from the raw FASTQ. This is the parent ID for split fragments. |
| `raw_fastq_file` | Raw FASTQ file containing the original read. |
| `raw_read_length` | Length of the original unfiltered read. |
| `forward_patternName` | Name of the best seqkit-detected forward primer/UMI hit on the raw read. Blank when no forward hit was found. |
| `forward_strand` | Strand of the best forward hit. |
| `forward_start` | Start coordinate of the best forward hit. |
| `forward_end` | End coordinate of the best forward hit. |
| `forward_mismatch_count` | Number of mismatches in the best forward hit. |
| `forward_primer_identity` | Best forward-hit identity score. |
| `reverse_patternName` | Name of the best seqkit-detected reverse primer/UMI hit on the raw read. Blank when no reverse hit was found. |
| `reverse_strand` | Strand of the best reverse hit. |
| `reverse_start` | Start coordinate of the best reverse hit. |
| `reverse_end` | End coordinate of the best reverse hit. |
| `reverse_mismatch_count` | Number of mismatches in the best reverse hit. |
| `reverse_primer_identity` | Best reverse-hit identity score. |
| `entered_amplicon_sorter` | `yes` when an original or split processed read was present in the amplicon_sorter input map; otherwise `no`. |
| `processed_read_id` | ID passed to amplicon_sorter. Split reads are `split-by-adapter-left-<UUID>` or `split-by-adapter-right-<UUID>`. |
| `processed_fragment` | `original`, `left`, or `right`. Blank when the raw read did not enter amplicon_sorter. |
| `processed_read_length` | Length after filtering, trimming, and possible adapter splitting. |
| `amplicon_sorter_group` | Exact amplicon_sorter subgroup, for example `1_0`. Blank for ungrouped reads. |
| `amplicon_sorter_group_consensus_id` | Consensus FASTA/Kraken record ID corresponding to the subgroup. |
| `group_consensus_kraken_taxid` | Kraken taxid assigned to the subgroup consensus and propagated to all reads in that subgroup. |
Reduced read-tracking table columns:

| Column | Meaning |
|---|---|
| `raw_ont_read_id` | Original raw ONT read ID. |
| `forward_patternName` | Best seqkit-detected forward primer/UMI hit. |
| `reverse_patternName` | Best seqkit-detected reverse primer/UMI hit. |
| `group_consensus_kraken_taxid` | Kraken taxid propagated from the amplicon_sorter group consensus. |

Taxid `0` means the consensus record was sent to Kraken2 and classified as unclassified. A blank taxid means no corresponding group consensus Kraken record was available.

## 10. Evaluate tag tracking and taxid results

Notebook:

```text
code/evaluate_tagtracking_taxid_reduced.ipynb
```

The notebook is parameterized by the `SELECTED_BARCODE` environment variable. It reads:

```text
data/mergedseq5/barcodeXX/barcodeXX_read_primer_group_taxid_reduced.tsv
data/mergedseq5/barcodeXX/barcodeXX_insecta_read_primer_group_taxid_reduced.tsv
```

Run one barcode:

```bash
mkdir -p code/executed_notebooks logs

SELECTED_BARCODE=barcode01 jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir code/executed_notebooks \
  --output evaluate_tagtracking_taxid_reduced_barcode01.ipynb \
  --ExecutePreprocessor.timeout=-1 \
  code/evaluate_tagtracking_taxid_reduced.ipynb \
  2>&1 | tee logs/evaluate_tagtracking_barcode01.log
```

Run all barcodes:

```bash
mkdir -p code/executed_notebooks logs

for n in $(seq -w 1 24); do
  barcode="barcode${n}"
  echo "Running ${barcode}"

  SELECTED_BARCODE="${barcode}" jupyter nbconvert \
    --execute \
    --to notebook \
    --output-dir code/executed_notebooks \
    --output "evaluate_tagtracking_taxid_reduced_${barcode}.ipynb" \
    --ExecutePreprocessor.timeout=-1 \
    code/evaluate_tagtracking_taxid_reduced.ipynb \
    > "logs/evaluate_tagtracking_${barcode}.log" 2>&1

  if [ $? -ne 0 ]; then
    echo "FAILED: ${barcode}"
    tail -80 "logs/evaluate_tagtracking_${barcode}.log"
    break
  fi
done
```

Monitor running notebooks:

```bash
ps -ef | grep -E 'jupyter|nbconvert|ipykernel' | grep -v grep
```

Check completed executed notebook copies:

```bash
ls -lh code/executed_notebooks/evaluate_tagtracking_taxid_reduced_barcode*.ipynb
```

Check recent exported tables:

```bash
find data/mergedseq5 -name '*target_species_detection.tsv' -printf '%TY-%Tm-%Td %TH:%TM %p\n' | sort | tail
```

Evaluation outputs are written beside each barcode's reduced tracking input:

```text
data/mergedseq5/barcodeXX/
  barcodeXX_top10_taxa_per_forward_tag.tsv
  barcodeXX_unsorted_taxa.tsv
  barcodeXX_primer_detection_classification_summary.tsv
  barcodeXX_target_species_detection.tsv
  barcodeXX_taxids_across_50pct_forward_primers.tsv
  barcodeXX_target_species_taxids_below_50pct_forward_primers.tsv
  barcodeXX_expected_species_not_detected.tsv
  figures/
    barcodeXX_primer_detection_classification.png
    barcodeXX_reads_per_forward_primer.png
    barcodeXX_reads_per_reverse_primer.png
    barcodeXX_unique_taxids_per_forward_primer.png
    barcodeXX_insecta_reads_vs_taxid_diversity.png
    barcodeXX_top10_taxids_per_forward_primer.png
```

The executed notebook copies go to:

```text
code/executed_notebooks/
```

The detailed column documentation for these exported tables is in:

```text
data/mergedseq5/README.md
```

## 11. Generate Sankey diagrams from target-species detections

Notebook:

```text
code/sankey_target_species_detection.ipynb
```

Run the notebook from the repository root after
`evaluate_tagtracking_taxid_reduced.ipynb` has created the per-barcode
`barcodeXX_target_species_detection.tsv` files:

```bash
mkdir -p code/executed_notebooks logs

DATE_TAG="$(date +%Y%m%d)"
TIME_TAG="$(date +%Y%m%d_%H%M%S)"

conda run -n saproxdna-evaluate jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir code/executed_notebooks \
  --output "sankey_target_species_detection_mergedseq5_${DATE_TAG}.ipynb" \
  --ExecutePreprocessor.timeout=-1 \
  code/sankey_target_species_detection.ipynb \
  > "logs/sankey_target_species_detection_${TIME_TAG}.log" 2>&1
```

The command executes every cell non-interactively, saves an executed notebook
under `code/executed_notebooks/`, and writes the run log under `logs/`. Check
the log if the command exits unsuccessfully:

```bash
tail -100 "$(find logs -maxdepth 1 -type f \
  -name 'sankey_target_species_detection_*.log' | sort | tail -n 1)"
```

By default, the notebook processes every matching barcode table below
`data/mergedseq5/`. To process only selected barcodes, set
`SELECTED_BARCODES` in the notebook settings cell before running it, for
example:

```python
SELECTED_BARCODES = ['barcode01', 'barcode06']
```

The full downstream wrapper runs the same notebook automatically as its final
step:

```bash
bash code/run_full_downstream_pipeline.sh
```

The wrapper command reruns the preceding downstream steps as well; use the
standalone `nbconvert` command above when the target-species detection tables
already exist and only the Sankey outputs are needed.

This notebook creates one Sankey/alluvial figure per barcode from:

```text
data/mergedseq5/barcodeXX/barcodeXX_target_species_detection.tsv
```

It writes an expanded trust-annotated table:

```text
data/mergedseq5/barcodeXX/barcodeXX_target_species_detection_trusted.tsv
```

The final notebook section also combines trusted rows across selected barcodes into:

```text
data/mergedseq5/mergedseq5_target_species_detection_trusted.tsv
```

It also writes two Sankey figures:

```text
data/mergedseq5/barcodeXX/figures/barcodeXX_target_detection_sankey.png
data/mergedseq5/barcodeXX/figures/barcodeXX_target_detection_sankey_trusted.png
```

The same notebook also reads `barcodeXX_target_species_taxids_below_50pct_forward_primers.tsv` and applies the same Sankey and trust-coloring rules. For that filtered table, it does not export a new expanded/trusted TSV; it only writes:

```text
data/mergedseq5/barcodeXX/figures/barcodeXX_target_species_taxids_below_50pct_forward_primers_sankey.png
data/mergedseq5/barcodeXX/figures/barcodeXX_target_species_taxids_below_50pct_forward_primers_sankey_trusted.png
```

The plotted flow is:

```text
expected_forward_primer -> detected_forward_primer -> detected_species -> target_species
```

Rules:

- The notebook uses the full `barcodeXX_target_species_detection.tsv` table, not the below-50% filtered table.
- `expected_forward_primer` is taken directly from the row's `expected_forward` column. Missing or empty values are grouped as `none`. Multiple semicolon-separated expected primers are split.
- `target_species` is taken directly from the same row's `target_species` column. It is not globally re-derived from the expected forward primer. Missing or empty values are grouped as `none`.
- `detected_species` is taken from `taxid_detected_name`. Missing or empty values are grouped as `none`.
- `detected_forward_primer` is taken from the seqkit-detected `forward_primer` column. Missing or empty detected forward primers are grouped as `none`.
- `target_sxdcodes` and `expected_forward` are expanded together. If both contain the same number of semicolon-separated values, they are paired by position so the first SxD code remains linked to the first expected forward primer, and so on. If one side has one value and the other has multiple, the single value is repeated. Otherwise, all combinations are used and the row is flagged in `expansion_note`.
- If one row expands to multiple expected-primer/detected-primer/species combinations, the row's read count is split evenly across those combinations as `sankey_read_count` so the total flow remains conserved.
- Link weights are based on read-count columns in `barcodeXX_target_species_detection.tsv`: `taxid_primer_read_count` for detected-forward-primer rows and `no_forward_primer_read_count` for no-forward-primer rows, with fallback to `taxid_read_count` and then `target_species_read_count`.
- Detected species and target species blocks are sorted alphabetically, with `none` at the end.
- Expected and detected forward-primer blocks are initialized alphanumerically/naturally, then refined with weighted barycentric sweeps to reduce path overlap/crossings. This is a layout-only operation and does not change link weights.
- The diagram is exploratory. It shows expected-primer, expected-target, detected-taxid, and detected-primer associations but does not prove that a detected taxid is the exact expected species.

Trust rules used for `barcodeXX_target_species_detection_trusted.tsv` and the green trusted Sankey:

- After expansion, primer-detected rows are trusted only when `expected_forward_primer` and `detected_forward_primer` are the same.
- If a `target_species + expected_forward_primer` pair has multiple detected forward primers, only rows where the detected primer matches the expected primer are trusted; mismatching detected primers are not trusted.
- If a target species has multiple expected forward primers, each expanded expected/detected primer row is evaluated independently.
- Rows with no detected forward primer are trusted only when `target_genus` and `taxid_detected_genus` match and the detected species is unique to one `target_species + expected_forward_primer` pair in that barcode.
- Rows with no detected forward primer are not trusted when the detected species appears in multiple target/species expected-primer pairs.
- The trusted Sankey colors trusted links in green and untrusted links in grey.

### Build annotated noNUMT consensus FASTAs for phylogenetics

After the Sankey notebook has written
`barcodeXX_target_species_detection_trusted.tsv` for every selected barcode,
run:

```bash
conda run -n saproxdna-evaluate python code/build_phylo_consensus_fastas.py \
  --merged-dir data/mergedseq5 \
  --barcode-start 1 \
  --barcode-end 24
```

The script reads, for each barcode:

```text
barcodeXX/filtered/combined/*consensussequences_noNUMT_*filter.fasta
barcodeXX/filtered/combined/*consensussequences_noNUMT_*filter.fasta.kraken
barcodeXX/barcodeXX_target_species_detection_trusted.tsv
```

It uses the Kraken database scientific-name table at
`/mnt/d/kraken2DB/k2_core_nt_2025/names.dmp` by default. Supply a different
file with `--names-dmp` if another Kraken database was used.

Every noNUMT consensus is retained and receives one of these annotations:

- `Kraken_scientific_name_trusted` when the consensus is classified and its
  taxID has at least one `trusted=True` association in that barcode's trusted
  table.
- `Kraken_scientific_name_notrusted` when Kraken classified the consensus but
  its taxID has no trusted association in that barcode.
- `noid` when the consensus is absent from the Kraken output, Kraken reports it
  as unclassified, or its taxID is `0`.

The original consensus header is sanitized by replacing each run of
non-alphanumeric characters, including spaces, with `_`. The barcode name is
also appended before the Kraken name and trust status, for example:

```text
>consensus_asinput_0_0_1545_barcode01_Loewia_foeda_trusted
```

Per-barcode outputs:

```text
data/mergedseq5/barcodeXX/phylo/barcodeXX_consensus_noNUMT_phylo.fasta
```

The script validates that sanitized headers are unique, then concatenates the
per-barcode FASTAs into:

```text
data/mergedseq5/consensus/mergedseq5_consensus_noNUMT_phylo.fasta
```

Outputs are written atomically, so a failed or interrupted run does not replace
a previously completed FASTA. `code/run_full_downstream_pipeline.sh` runs this
export automatically after the Sankey notebook.

## 12. Recommended full run order

For a complete fresh run after raw reads have been placed in barcode folders:

1. Run `code/run_metabarcoding_pipeline.sh` for all barcode folders. This runs fastplong, amplicon_sorter, noNUMT consensus filtering, Kraken2, and the Kraken read map.
2. Run `code/length_check_amplicon_sorter.sh` after amplicon_sorter and Kraken2/read-map generation have produced `filtered/combined/asinput_*_*.fasta`, `asinput_consensussequences_noNUMT_200filter.read_map.tsv`, `asinput_consensussequences_noNUMT_200filter.fasta.kraken`, and `asinput_consensussequences_noNUMT_200filter.fasta.report` files.
3. Run `code/run_primer_umi_tracking.sh` for the parent merged-sequence folder.
4. Run `code/primer_umi_seqkit_parser.ipynb` to create combined primer/UMI best-hit tables under `data/primertrack/`.
5. Run `code/primer_diversity_amplicon_sorter.sh` after the seqkit best-hit table exists.
6. Run `code/kraken_parser.ipynb` to create target-clade Kraken parsed tables under `data/kraken_parsed/`.
7. Run `code/build_barcode_read_tracking_table.py` for each barcode.
8. Run `code/evaluate_tagtracking_taxid_reduced.ipynb` for each barcode using `SELECTED_BARCODE`.
9. Run `code/sankey_target_species_detection.ipynb` if Sankey diagrams are needed from the exported target-species detection tables.
10. Run `code/build_phylo_consensus_fastas.py` to annotate the noNUMT consensus sequences and create the per-barcode and combined phylogenetic FASTAs.

## 13. Rerun guidance

Use `FORCE=1` only when you intentionally want to redo completed fastplong/amplicon_sorter/Kraken outputs.

Use `FORCE_KRAKEN=1` when the Kraken database, Kraken parameters, or read-map generation changed but filtering and amplicon_sorter outputs can be kept.

If primer FASTA or primer mismatch settings change, rerun:

1. `code/run_primer_umi_tracking.sh`
2. `code/primer_umi_seqkit_parser.ipynb`
3. `code/build_barcode_read_tracking_table.py`
4. `code/evaluate_tagtracking_taxid_reduced.ipynb`

If taxonomy parsing or target clade selection changes, rerun:

1. `code/kraken_parser.ipynb`
2. `code/build_barcode_read_tracking_table.py`
3. `code/evaluate_tagtracking_taxid_reduced.ipynb`

If only plotting/table interpretation changed in the evaluation notebook, rerun only:

```text
code/evaluate_tagtracking_taxid_reduced.ipynb
```

## 14. Temporary files

amplicon_sorter uses:

```text
/tmp/metabarcoding_amplicon_sorter/
```

Successful runs copy usable output back to the barcode folder. Failed runs may leave temporary folders for inspection. Do not remove temporary folders while amplicon_sorter is running.

List temporary amplicon_sorter folders:

```bash
find /tmp/metabarcoding_amplicon_sorter -mindepth 1 -maxdepth 1 -type d -print
```

Remove only failed-run directories that you no longer need for inspection.

Do not remove filtered reads, amplicon_sorter outputs, Kraken2 outputs, or seqkit outputs unless you intentionally want to rerun those steps.

Common cleanup candidates are:

```text
code/__pycache__/
tmp/
```

Older parser outputs without the current dated naming scheme may also be stale, for example undated files in `data/kraken_parsed/`. Confirm that the dated files you need exist before deleting older parser summaries.
