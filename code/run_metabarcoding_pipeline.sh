#!/usr/bin/env bash
# Module 1 per sample: primers and tags (cutadapt + primer_tags.py) -> quality filter (fastplong)
# -> cluster each specimen at 97% (cluster_specimens.py) -> classify the clusters (VSEARCH SINTAX).
set -euo pipefail
CONFIG_FILE="${CONFIG_FILE:-config/pipeline.env}"
[[ -f "$CONFIG_FILE" ]] || { echo "Missing settings file: $CONFIG_FILE (run from the project root)" >&2; exit 1; }
source "$CONFIG_FILE"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"

prepare_database() (
    # Called only when classification is needed, before starting sample workers.
    [[ ! -s "$SINTAX_DB" ]] || return 0
    mkdir -p "$(dirname "$SINTAX_DB")"
    local temp_ref="" temp_db=""
    trap '[[ -z "$temp_ref" ]] || rm -f -- "$temp_ref"; [[ -z "$temp_db" ]] || rm -f -- "$temp_db"' EXIT
    if [[ ! -s "$SINTAX_REFERENCE" ]]; then
        [[ -s "${SINTAX_REFERENCE}.gz" ]] || { echo "Missing reference: $SINTAX_REFERENCE[.gz]" >&2; exit 1; }
        temp_ref="$(mktemp "${SINTAX_REFERENCE}.XXXXXX")"
        gzip -dc -- "${SINTAX_REFERENCE}.gz" > "$temp_ref"
        mv -- "$temp_ref" "$SINTAX_REFERENCE"; temp_ref=""
    fi
    temp_db="$(mktemp "${SINTAX_DB}.XXXXXX")"
    vsearch --makeudb_usearch "$SINTAX_REFERENCE" --output "$temp_db"
    mv -- "$temp_db" "$SINTAX_DB"; temp_db="" )

cleaned_ready() {
    [[ -s "$1/filtered/asinput.fastq.gz" && -s "$1/filtered/primer_tags.tsv.gz" ]]; }

clean_sample() (
    # Raw reads -> primer/tag table + primer-trimmed, quality-filtered reads.
    local folder="$1" filtered="$1/filtered" sample work fifo_dir pid
    sample="$(basename "$folder")"
    if [[ -e "$filtered" ]] && ! cleaned_ready "$folder"; then
        [[ "$FORCE" == 1 ]] || { echo "[$sample] $filtered is incomplete or was made by an older pipeline; set FORCE=1 to replace it" >&2; exit 1; }
    fi
    rm -rf -- "$filtered"
    mkdir -p "$filtered/logs"
    shopt -s nullglob
    local files=("$folder"/*.fastq.gz "$folder"/*.fq.gz)
    shopt -u nullglob
    [[ "${#files[@]}" -gt 0 ]] || { echo "[$sample] no raw *.fastq.gz files directly inside the folder" >&2; exit 1; }
    work="$(mktemp -d "$filtered/.work.XXXXXX")"
    # Windows drives mounted in WSL (/mnt/c) cannot hold a named pipe, so it lives in the Linux temp folder.
    fifo_dir="$(mktemp -d "${TMPDIR:-/tmp}/primer_tags.XXXXXX")"
    trap 'rm -rf -- "$work" "$fifo_dir"' EXIT
    mkfifo "$fifo_dir/info.fifo"
    echo "[$sample] cutadapt + tag scoring on ${#files[@]} raw FASTQ files"
    "$PYTHON" "$CODE_DIR/primer_tags.py" --primers "$PRIMER_FASTA" --info "$fifo_dir/info.fifo" \
        --output "$work/primer_tags.tsv.gz" --max-tag-errors "$TAG_MAX_ERRORS" \
        --min-margin "$TAG_MIN_MARGIN" --min-observed "$TAG_MIN_OBSERVED" 2> "$filtered/logs/primer_tags.log" &
    pid=$!
    # Multi-core cutadapt needs an uncompressed stream; gzip members concatenate safely.
    # Reads shorter than fastplong's minimum are dropped here: fastplong would remove them
    # anyway, and it stops reading its input at the first empty read (nothing between the
    # primers). The tag table still lists every raw read.
    if ! gzip -dc -- "${files[@]}" | "$CUTADAPT_BIN" -j "$PRIMER_THREADS" --revcomp -g "$PRIMER_ADAPTER" \
            --minimum-length "$FASTPLONG_LENGTH_REQUIRED" \
            --info-file "$fifo_dir/info.fifo" -o "$work/primer_trimmed.fastq.gz" - > "$filtered/logs/cutadapt.log" 2>&1; then
        kill "$pid" 2>/dev/null || true
        echo "[$sample] cutadapt failed; see $filtered/logs/cutadapt.log" >&2; exit 1
    fi
    wait "$pid" || { echo "[$sample] tag scoring failed; see $filtered/logs/primer_tags.log" >&2; exit 1; }
    cat "$filtered/logs/primer_tags.log"
    local fastplong_cmd=(fastplong -i "$work/primer_trimmed.fastq.gz" -o "$work/asinput.fastq.gz"
        --unqualified_percent_limit "$FASTPLONG_UNQUALIFIED_PERCENT_LIMIT"
        --qualified_quality_phred "$FASTPLONG_QUALIFIED_QUALITY_PHRED"
        --mean_qual "$FASTPLONG_MEAN_QUAL" --length_required "$FASTPLONG_LENGTH_REQUIRED"
        --length_limit "$FASTPLONG_LENGTH_LIMIT" --trimming_extension "$FASTPLONG_TRIMMING_EXTENSION"
        --adapter_fasta "$ADAPTER_FASTA" --start_adapter "$FASTPLONG_START_ADAPTER"
        --end_adapter "$FASTPLONG_END_ADAPTER" --thread "$FASTPLONG_THREADS"
        -j "$filtered/fastplong.json" -h "$filtered/fastplong.html")
    echo "[$sample] fastplong quality filtering"
    { printf 'Command:'; printf ' %q' "${fastplong_cmd[@]}"; printf '\n\n'; } > "$filtered/logs/fastplong.log"
    # fastplong prints some trimming details only to a terminal; script records them.
    script -q -e -a -c "$(printf '%q ' "${fastplong_cmd[@]}")" "$filtered/logs/fastplong.log" > /dev/null
    # Guard against fastplong silently reading only part of its input.
    "$PYTHON" - "$filtered/logs/cutadapt.log" "$filtered/fastplong.json" <<'PY' || exit 1
import json, re, sys
written = int(re.search(r'Reads written \(passing filters\): +([0-9,]+)', open(sys.argv[1]).read())[1].replace(',', ''))
read = json.load(open(sys.argv[2]))['summary']['before_filtering']['total_reads']
if read != written:
    sys.exit(f'fastplong read {read} of the {written} reads written by cutadapt; see {sys.argv[2]}')
PY
    mv -- "$work/primer_tags.tsv.gz" "$filtered/primer_tags.tsv.gz"
    mv -- "$work/asinput.fastq.gz" "$filtered/asinput.fastq.gz" )

clusters_ready() {
    # Clusters exist and were made from the current cleaned reads.
    [[ -s "$1/filtered/clusters.fasta" && -s "$1/filtered/cluster_members.tsv.gz" &&
       ! "$1/filtered/asinput.fastq.gz" -nt "$1/filtered/clusters.fasta" ]]; }

classification_output() {
    local folder="$1" sample dataset
    sample="$(basename "$folder")"
    dataset="$(basename "$(dirname "$folder")")"
    printf '%s/%s/%s_%s.tsv\n' "${SINTAX_OUTPUT_ROOT:-data/sintax/$dataset}" "$sample" "$sample" "$SINTAX_SUFFIX"; }

needs_cleaning() { [[ "$FORCE" == 1 ]] || ! cleaned_ready "$1"; }

needs_clustering() { needs_cleaning "$1" || ! clusters_ready "$1"; }

needs_sintax() {
    local output
    output="$(classification_output "$1")"
    # Reclassify when the clusters are newer than the existing classification.
    needs_clustering "$1" || [[ "$FORCE_SINTAX" == 1 || ! -f "$output" || "$1/filtered/clusters.fasta" -nt "$output" ]]; }

run_folder() (
    local folder="$1" sample output tmp=""
    sample="$(basename "$folder")"
    output="$(classification_output "$folder")"
    if needs_cleaning "$folder"; then clean_sample "$folder"; else echo "[$sample] reuse cleaned reads"; fi
    if needs_clustering "$folder"; then
        echo "[$sample] cluster each specimen at $CLUSTER_IDENTITY identity"
        "$PYTHON" "$CODE_DIR/cluster_specimens.py" "$folder" --identity "$CLUSTER_IDENTITY" --threads "$CLUSTER_THREADS" |
            tee "$folder/filtered/logs/clustering.log"
    else
        echo "[$sample] reuse clusters"
    fi
    [[ "$CLEAN_ONLY" != 1 ]] || return 0
    if ! needs_sintax "$folder"; then echo "[$sample] reuse $output"; return 0; fi
    mkdir -p "$(dirname "$output")"
    tmp="$(mktemp "${output}.XXXXXX")"
    trap '[[ -z "$tmp" ]] || rm -f -- "$tmp"' EXIT
    vsearch --sintax "$folder/filtered/clusters.fasta" --db "$SINTAX_DB" \
        --sintax_random --strand both --sintax_cutoff "$SINTAX_CUTOFF" \
        --threads "$SINTAX_THREADS" --randseed "$SINTAX_RANDSEED" \
        --maxseqlength "$SINTAX_MAX_LENGTH" --minseqlength "$SINTAX_MIN_LENGTH" \
        --tabbedout "$tmp"
    mv -- "$tmp" "$output"; tmp="" )

main() {
    if [[ "$#" -eq 0 || "${1:-}" == --help ]]; then
        echo 'Usage: bash code/run_metabarcoding_pipeline.sh [--dry-run] FOLDER [FOLDER ...]'
        echo 'Each FOLDER is a dataset parent or a barcodeNN/unclassified sample. Bare dataset names resolve under data/.'
        return 0
    fi
    local dry_run=0
    if [[ "$1" == --dry-run ]]; then dry_run=1; shift; fi
    [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "MAX_JOBS must be a positive integer" >&2; return 2; }
    local folder needs_db=0 needs_clean=0 failed=0 pid discovered
    local -a folders=() pids=()
    discovered="$("$PYTHON" "$CODE_DIR/sample_folders.py" "$@")" || return 1
    mapfile -t folders <<< "$discovered"
    for folder in "${folders[@]}"; do
        if needs_cleaning "$folder"; then needs_clean=1; echo "[$(basename "$folder")] clean raw reads -> filtered/"; fi
        if needs_clustering "$folder"; then echo "[$(basename "$folder")] cluster specimens -> filtered/clusters.fasta"; fi
        if [[ "$CLEAN_ONLY" == 1 ]]; then continue; fi
        if needs_sintax "$folder"; then needs_db=1; echo "[$(basename "$folder")] classify -> $(classification_output "$folder")"
        else echo "[$(basename "$folder")] reuse $(classification_output "$folder")"; fi
    done
    [[ "$dry_run" != 1 ]] || return 0
    if [[ "$needs_clean" == 1 ]]; then
        for file in "$PRIMER_FASTA" "$ADAPTER_FASTA"; do [[ -f "$file" ]] || { echo "Missing $file" >&2; return 1; }; done
        command -v fastplong >/dev/null || { echo "fastplong not found; activate the taxonomy environment" >&2; return 1; }
        CUTADAPT_BIN="${CUTADAPT_BIN:-$(command -v cutadapt || "${CONDA_EXE:-conda}" run -n "$CUTADAPT_ENV" which cutadapt)}" ||
            { echo "cutadapt not found; create the '$CUTADAPT_ENV' environment (docs/setup.md)" >&2; return 1; }
        PRIMER_ADAPTER="$("$PYTHON" "$CODE_DIR/primer_tags.py" --primers "$PRIMER_FASTA" --print-adapter "$PRIMER_CORE_MAX_ERRORS")"
        export CUTADAPT_BIN PRIMER_ADAPTER
    fi
    command -v vsearch >/dev/null || { echo "vsearch not found; activate the taxonomy environment" >&2; return 1; }
    if [[ "$needs_db" == 1 ]]; then prepare_database; fi
    for folder in "${folders[@]}"; do
        run_folder "$folder" &
        pids+=("$!")
        if [[ "${#pids[@]}" -ge "$MAX_JOBS" ]]; then
            for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
            pids=()
        fi
    done
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    return "$failed"; }
main "$@"
