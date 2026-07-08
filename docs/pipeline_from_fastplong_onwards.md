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
| `FASTPLONG_LENGTH_REQUIRED` | fastplong minimum read length. |
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
| `PRIMER_UMI_MAX_MISMATCH` | Maximum mismatches allowed by `seqkit locate`. |
| `PRIMER_UMI_SEQKIT_THREADS` | Threads used by seqkit primer tracking. |

Activate the environment that has fastplong, amplicon_sorter dependencies, Kraken2, seqkit, Jupyter, pandas, and matplotlib available. In the current work this has been the `saproxdna` environment for the shell pipeline and `saproxdna-evaluate` for notebook evaluation.

## 2. Run fastplong, amplicon_sorter, Kraken2, and Kraken read-map generation

Script:

```text
code/run_metabarcoding_pipeline.sh
```

This script runs, per barcode:

1. fastplong quality/filter/adapter trimming.
2. Combine filtered FASTQ files into `filtered/asinput.fastq.gz`.
3. amplicon_sorter using `code/amplicon_sorter_mod.py`.
4. Build `kraken2_combined_input.fasta` from amplicon_sorter consensus and unique-read FASTA files.
5. Build `kraken2_combined_input.read_map.tsv`, linking Kraken record IDs back to amplicon_sorter read numbers and ONT read IDs.
6. Kraken2 classification of the combined FASTA.

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
    asinput_nogroup_unique.fasta
    kraken2_combined_input.fasta
    kraken2_combined_input.read_map.tsv
    kraken2_combined_input.fasta.kraken
    kraken2_combined_input.fasta.report
```

### Important read-map outputs

`filtered/combined/amplicon_sorter_read_map.tsv` links amplicon_sorter internal read numbers to ONT read IDs.

| Column | Meaning |
|---|---|
| `amplicon_sorter_read_number` | Numeric read identifier assigned internally by amplicon_sorter. |
| `ont_read_id` | Original or processed ONT read ID. Split reads can be named `split-by-adapter-left-<raw_id>` or `split-by-adapter-right-<raw_id>`. |
| `read_length` | Length of the sequence passed into amplicon_sorter. |

`filtered/combined/kraken2_combined_input.read_map.tsv` links every Kraken input record to amplicon_sorter context.

| Column | Meaning |
|---|---|
| `merge_folder` | Parent sequencing merge/run folder, for example `mergedseq5`. |
| `barcode` | Barcode folder, for example `barcode01`. |
| `kraken_read_id` | FASTA record ID sent to Kraken2. Numeric IDs are individual reads; consensus IDs start with `consensus`. |
| `record_type` | `read`, `consensus`, or `other`. |
| `amplicon_sorter_file` | amplicon_sorter FASTA file containing this Kraken input sequence. |
| `amplicon_sorter_group` | Group/file stem associated with this record, for example `asinput_1_unique`, `asinput_nogroup_unique`, or `asinput_consensussequences`. |
| `amplicon_sorter_read_number` | amplicon_sorter internal read number for individual-read records. Blank for consensus records. |
| `ont_read_id` | ONT read ID for individual-read records. Blank for consensus records. |
| `ont_read_length` | ONT/processed sequence length for individual-read records. |
| `kraken_sequence_length` | Length of the exact sequence sent to Kraken2. |

## 3. Run primer/UMI tracking with seqkit

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

## 4. Parse primer/UMI hits into best-hit tables

Notebook:

```text
code/primer_umi_seqkit_parser.ipynb
```

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

## 5. Parse Kraken2 outputs

Notebook:

```text
code/kraken_parser.ipynb
```

This notebook reads Kraken2 outputs and `kraken2_combined_input.read_map.tsv` files, attaches barcode/run metadata, and exports parsed Kraken/taxonomy tables under:

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

The read-tracking script currently uses the Insecta reads file:

```text
data/kraken_parsed/mergedseq5_insecta_kraken_reads_20260630.csv
```

Required columns in that file:

| Column | Meaning |
|---|---|
| `readID` | Kraken input record ID. This can be a consensus ID or individual-read ID. |
| `merge_folder` | Merge/run folder. |
| `barcode` | Barcode folder. |

Additional taxonomy columns are produced by the parser and used for reporting, but `build_barcode_read_tracking_table.py` only needs the three columns above to intersect reads/groups with a target clade such as Insecta.

## 6. Build per-read primer/group/taxid tracking tables

Script:

```text
code/build_barcode_read_tracking_table.py
```

This script combines:

- raw ONT read IDs and raw read lengths;
- seqkit best forward/reverse primer hits;
- amplicon_sorter processed read IDs and split-read fragments;
- amplicon_sorter subgroup membership;
- Kraken2 group-consensus taxids;
- Kraken2 individual-read taxids;
- the Insecta subset from `kraken_parser.ipynb`.

Run one barcode:

```bash
python code/build_barcode_read_tracking_table.py \
  --barcode-dir data/mergedseq5/barcode01 \
  --primer-best-hits data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_20260630.txt \
  --insecta-kraken-reads data/kraken_parsed/mergedseq5_insecta_kraken_reads_20260630.csv \
  --output data/mergedseq5/barcode01/barcode01_read_primer_group_taxid.tsv
```

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
| `individual_read_kraken_input_file` | amplicon_sorter FASTA used as individual-read Kraken input, for example `asinput_nogroup_unique.fasta`, `asinput_0_unique.fasta`, or `asinput_1_unique.fasta`. |
| `individual_read_kraken_taxid` | Kraken taxid assigned directly to the individual unique-read sequence. |

Reduced read-tracking table columns:

| Column | Meaning |
|---|---|
| `raw_ont_read_id` | Original raw ONT read ID. |
| `forward_patternName` | Best seqkit-detected forward primer/UMI hit. |
| `reverse_patternName` | Best seqkit-detected reverse primer/UMI hit. |
| `group_consensus_kraken_taxid` | Kraken taxid propagated from the amplicon_sorter group consensus. |
| `individual_read_kraken_taxid` | Kraken taxid assigned directly to an individual ungrouped/unique read. |

Taxid `0` means the record was sent to Kraken2 and classified as unclassified. A blank taxid means no corresponding Kraken record was available.

## 7. Evaluate tag tracking and taxid results

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

## 8. Recommended full run order

For a complete fresh run after raw reads have been placed in barcode folders:

1. Run `code/run_metabarcoding_pipeline.sh` for all barcode folders. This runs fastplong, amplicon_sorter, Kraken2, and the Kraken read map.
2. Run `code/run_primer_umi_tracking.sh` for the parent merged-sequence folder.
3. Run `code/primer_umi_seqkit_parser.ipynb` to create combined primer/UMI best-hit tables under `data/primertrack/`.
4. Run `code/kraken_parser.ipynb` to create target-clade Kraken parsed tables under `data/kraken_parsed/`.
5. Run `code/build_barcode_read_tracking_table.py` for each barcode.
6. Run `code/evaluate_tagtracking_taxid_reduced.ipynb` for each barcode using `SELECTED_BARCODE`.

## 9. Rerun guidance

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

## 10. Temporary files

amplicon_sorter uses:

```text
/tmp/metabarcoding_amplicon_sorter/
```

Successful runs copy usable output back to the barcode folder. Failed runs may leave temporary folders for inspection. Do not remove temporary folders while amplicon_sorter is running.

List temporary amplicon_sorter folders:

```bash
find /tmp/metabarcoding_amplicon_sorter -mindepth 1 -maxdepth 1 -type d -print
```

