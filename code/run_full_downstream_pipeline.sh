#!/usr/bin/env bash
set -euo pipefail

SHELL_ENV="${SHELL_ENV:-saproxdna}"
SEQKIT_ENV="${SEQKIT_ENV:-umitracker}"
NOTEBOOK_ENV="${NOTEBOOK_ENV:-saproxdna-evaluate}"
KRAKEN_NAMES_DMP="${KRAKEN_NAMES_DMP:-/mnt/d/kraken2DB/k2_core_nt_2025/names.dmp}"
MERGED_DIR="${MERGED_DIR:-data/mergedseq5}"
MERGE_NAME="$(basename "$MERGED_DIR")"
BARCODE_START="${BARCODE_START:-1}"
BARCODE_END="${BARCODE_END:-24}"
DATE_TAG="$(date +%Y%m%d)"
TIME_TAG="$(date +%Y%m%d_%H%M%S)"

mkdir -p logs code/executed_notebooks

echo "Running downstream pipeline"
echo "SHELL_ENV=$SHELL_ENV"
echo "SEQKIT_ENV=$SEQKIT_ENV"
echo "NOTEBOOK_ENV=$NOTEBOOK_ENV"
echo "KRAKEN_NAMES_DMP=$KRAKEN_NAMES_DMP"
echo "MERGED_DIR=$MERGED_DIR"
echo "BARCODES=barcode$(printf '%02d' "$BARCODE_START")..barcode$(printf '%02d' "$BARCODE_END")"
echo "TIME_TAG=$TIME_TAG"

barcode_dirs=()
for n in $(seq -w "$BARCODE_START" "$BARCODE_END"); do
  barcode_dir="${MERGED_DIR}/barcode${n}"
  if [[ ! -d "$barcode_dir" ]]; then
    echo "Missing barcode directory: $barcode_dir" >&2
    exit 1
  fi
  barcode_dirs+=("$barcode_dir")
done

echo "Step 0: fastplong, amplicon_sorter, noNUMT filtering, and Kraken2"
conda run --no-capture-output -n "$SHELL_ENV" bash code/run_metabarcoding_pipeline.sh "${barcode_dirs[@]}" \
  > "logs/run_metabarcoding_pipeline_${TIME_TAG}.log" 2>&1

echo "Step 1: primer/UMI tracking"
conda run --no-capture-output -n "$SEQKIT_ENV" bash code/run_primer_umi_tracking.sh "$MERGED_DIR" \
  > "logs/run_primer_umi_tracking_${TIME_TAG}.log" 2>&1

echo "Step 2: primer/UMI parser notebook"
SEQKIT_INPUT_FOLDER="$MERGED_DIR" conda run --no-capture-output -n "$NOTEBOOK_ENV" jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir code/executed_notebooks \
  --output "primer_umi_seqkit_parser_mergedseq5_${DATE_TAG}.ipynb" \
  --ExecutePreprocessor.timeout=-1 \
  code/primer_umi_seqkit_parser.ipynb \
  > "logs/primer_umi_seqkit_parser_${TIME_TAG}.log" 2>&1

echo "Step 3: primer diversity by amplicon_sorter file"
PRIMER_BEST_HITS="$(
  find data/primertrack \
    -maxdepth 1 \
    -type f \
    -name "${MERGE_NAME}_seqkit_primer_umi_best_hits_*.txt" \
    -printf '%T@\t%p\n' \
    | sort -n \
    | tail -n 1 \
    | cut -f 2-
)"
if [[ -z "$PRIMER_BEST_HITS" || ! -s "$PRIMER_BEST_HITS" ]]; then
  echo "Could not find a non-empty primer best-hit table for ${MERGE_NAME}" >&2
  exit 1
fi
conda run --no-capture-output -n "$SHELL_ENV" bash code/primer_diversity_amplicon_sorter.sh \
  "$MERGED_DIR" "$PRIMER_BEST_HITS" \
  > "logs/primer_diversity_amplicon_sorter_${TIME_TAG}.log" 2>&1

echo "Step 4: Kraken parser notebook"
conda run --no-capture-output -n "$NOTEBOOK_ENV" jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir code/executed_notebooks \
  --output "kraken_parser_mergedseq5_${DATE_TAG}.ipynb" \
  --ExecutePreprocessor.timeout=-1 \
  code/kraken_parser.ipynb \
  > "logs/kraken_parser_${TIME_TAG}.log" 2>&1

INSECTA_KRAKEN_READS="$(ls -t data/kraken_parsed/*insecta*reads* | head -1)"

echo "Using primer best hits: $PRIMER_BEST_HITS"
echo "Using Insecta Kraken reads: $INSECTA_KRAKEN_READS"

echo "Step 5: build barcode read-tracking tables"
for n in $(seq -w "$BARCODE_START" "$BARCODE_END"); do
  barcode="barcode${n}"
  echo "  $barcode"

  conda run --no-capture-output -n "$NOTEBOOK_ENV" python code/build_barcode_read_tracking_table.py \
    --barcode-dir "${MERGED_DIR}/${barcode}" \
    --primer-best-hits "$PRIMER_BEST_HITS" \
    --insecta-kraken-reads "$INSECTA_KRAKEN_READS" \
    --output "${MERGED_DIR}/${barcode}/${barcode}_read_primer_group_taxid.tsv" \
    > "logs/build_read_tracking_${barcode}_${TIME_TAG}.log" 2>&1
done

echo "Step 6: evaluate tag tracking notebooks"
for n in $(seq -w "$BARCODE_START" "$BARCODE_END"); do
  barcode="barcode${n}"
  echo "  $barcode"

  if [[ -f "${MERGED_DIR}/${barcode}/filtered/no_consensus.done" ]]; then
    echo "    no consensus sequences; skipping taxonomy evaluation notebook"
    continue
  fi

  SELECTED_BARCODE="$barcode" conda run --no-capture-output -n "$NOTEBOOK_ENV" jupyter nbconvert \
    --execute \
    --to notebook \
    --output-dir code/executed_notebooks \
    --output "evaluate_tagtracking_taxid_reduced_${barcode}_${DATE_TAG}.ipynb" \
    --ExecutePreprocessor.timeout=-1 \
    code/evaluate_tagtracking_taxid_reduced.ipynb \
    > "logs/evaluate_tagtracking_${barcode}_${TIME_TAG}.log" 2>&1
done

echo "Step 7: Sankey target-species detection notebook"
conda run --no-capture-output -n "$NOTEBOOK_ENV" jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir code/executed_notebooks \
  --output "sankey_target_species_detection_mergedseq5_${DATE_TAG}.ipynb" \
  --ExecutePreprocessor.timeout=-1 \
  code/sankey_target_species_detection.ipynb \
  > "logs/sankey_target_species_detection_${TIME_TAG}.log" 2>&1

echo "Step 8: annotated noNUMT consensus FASTAs for phylogenetics"
conda run --no-capture-output -n "$NOTEBOOK_ENV" python code/build_phylo_consensus_fastas.py \
  --merged-dir "$MERGED_DIR" \
  --barcode-start "$BARCODE_START" \
  --barcode-end "$BARCODE_END" \
  --names-dmp "$KRAKEN_NAMES_DMP" \
  > "logs/build_phylo_consensus_fastas_${TIME_TAG}.log" 2>&1

echo "Done. Logs use TIME_TAG=$TIME_TAG"
