#!/usr/bin/env bash
set -euo pipefail

MERGED_DIR="${1:-data/mergedseq5}"
OUTPUT_FILE="${2:-${MERGED_DIR}/mergedseq5_amplicon_sorter_length_check.tsv}"

mkdir -p "$(dirname "$OUTPUT_FILE")"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

printf 'barcode\tasinput_file\tcluster_group\tcluster_number\tgroup_number\tsequence_count\taverage_length\tlength_standard_deviation\tkraken_record_type\tkraken_record_count\tkraken_classified_count\tkraken_unclassified_count\tkraken_taxids\tkraken_taxon_names\tkraken_taxid_counts\tfile_path\n' > "$OUTPUT_FILE"

for barcode_dir in "$MERGED_DIR"/barcode*; do
    [[ -d "$barcode_dir/filtered/combined" ]] || continue

    barcode="$(basename "$barcode_dir")"
    combined_dir="$barcode_dir/filtered/combined"
    kraken_file="$combined_dir/kraken2_combined_input.fasta.kraken"
    kraken_map="$combined_dir/kraken2_combined_input.read_map.tsv"
    kraken_report="$combined_dir/kraken2_combined_input.fasta.report"
    kraken_summary="$TMP_DIR/${barcode}_kraken_by_asinput.tsv"

    if [[ -f "$kraken_file" && -f "$kraken_map" && -f "$kraken_report" ]]; then
        awk -F '\t' -v kraken_report="$kraken_report" -v kraken_file="$kraken_file" '
            BEGIN {
                taxid_name["0"] = "unclassified"
                taxid_name["NA"] = "NA"
            }
            FILENAME == kraken_report {
                taxid = $5
                name = $6
                gsub(/^[[:space:]]+/, "", name)
                gsub(/[[:space:]]+$/, "", name)
                if (taxid != "" && name != "") {
                    taxid_name[taxid] = name
                }
                next
            }
            FILENAME == kraken_file {
                kraken_status[$2] = $1
                kraken_taxid[$2] = $3
                next
            }
            FNR == 1 {
                next
            }
            {
                kraken_read_id = $3
                record_type = $4
                asinput_file = $5
                taxid = kraken_taxid[kraken_read_id]
                status = kraken_status[kraken_read_id]

                if (taxid == "") {
                    taxid = "NA"
                }
                if (status == "") {
                    status = "NA"
                }

                if (record_type == "consensus") {
                    group = kraken_read_id
                    sub(/^consensus_combinedrawreads_/, "", group)
                    sub(/\(.*/, "", group)
                    if (group != kraken_read_id && group != "") {
                        asinput_file = "asinput_" group ".fasta"
                    }
                    summary_record_type = "consensus"
                } else if (record_type == "read") {
                    summary_record_type = "unique_read"
                } else {
                    summary_record_type = record_type
                }

                record_count[asinput_file] += 1
                record_types[asinput_file, summary_record_type] = 1
                taxid_count[asinput_file, taxid] += 1
                if (taxid_order_seen[asinput_file, taxid] != 1) {
                    taxid_order_seen[asinput_file, taxid] = 1
                    taxid_order[asinput_file] = taxid_order[asinput_file] taxid ";"
                }
                if (status == "C" && taxid != "0" && taxid != "NA") {
                    classified_count[asinput_file] += 1
                } else {
                    unclassified_count[asinput_file] += 1
                }
            }
            END {
                for (asinput_file in record_count) {
                    record_type_list = ""
                    for (type_key in record_types) {
                        split(type_key, parts, SUBSEP)
                        if (parts[1] == asinput_file) {
                            record_type_list = record_type_list parts[2] ";"
                        }
                    }
                    sub(/;$/, "", record_type_list)

                    taxids = taxid_order[asinput_file]
                    sub(/;$/, "", taxids)
                    split(taxids, taxid_values, ";")
                    taxid_counts = ""
                    taxon_names = ""
                    for (taxid_index in taxid_values) {
                        taxid = taxid_values[taxid_index]
                        if (taxid != "") {
                            taxon_name = taxid_name[taxid]
                            if (taxon_name == "") {
                                taxon_name = "unknown_taxid_" taxid
                            }
                            taxid_counts = taxid_counts taxid ":" taxid_count[asinput_file, taxid] ";"
                            taxon_names = taxon_names taxon_name ";"
                        }
                    }
                    sub(/;$/, "", taxid_counts)
                    sub(/;$/, "", taxon_names)

                    printf "%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n",
                        asinput_file,
                        record_type_list,
                        record_count[asinput_file],
                        classified_count[asinput_file] + 0,
                        unclassified_count[asinput_file] + 0,
                        taxids,
                        taxon_names,
                        taxid_counts
                }
            }
        ' "$kraken_report" "$kraken_file" "$kraken_map" > "$kraken_summary"
    else
        : > "$kraken_summary"
    fi

    declare -A kraken_record_type_by_file=()
    declare -A kraken_record_count_by_file=()
    declare -A kraken_classified_count_by_file=()
    declare -A kraken_unclassified_count_by_file=()
    declare -A kraken_taxids_by_file=()
    declare -A kraken_taxon_names_by_file=()
    declare -A kraken_taxid_counts_by_file=()

    while IFS=$'\t' read -r summary_file record_type record_count classified_count unclassified_count taxids taxon_names taxid_counts; do
        kraken_record_type_by_file["$summary_file"]="$record_type"
        kraken_record_count_by_file["$summary_file"]="$record_count"
        kraken_classified_count_by_file["$summary_file"]="$classified_count"
        kraken_unclassified_count_by_file["$summary_file"]="$unclassified_count"
        kraken_taxids_by_file["$summary_file"]="$taxids"
        kraken_taxon_names_by_file["$summary_file"]="$taxon_names"
        kraken_taxid_counts_by_file["$summary_file"]="$taxid_counts"
    done < "$kraken_summary"

    find "$combined_dir" \
        -maxdepth 1 \
        -type f \
        -name 'asinput_*_*.fasta' \
        ! -iname '*alignment*fasta' \
        -print0 \
        | sort -z \
        | while IFS= read -r -d '' fasta_file; do
        asinput_file="$(basename "$fasta_file")"
        stem="${asinput_file%.fasta}"

        cluster_number='none'
        group_number='none'
        cluster_group="$stem"

        if [[ "$stem" =~ ^asinput_(.+)_unique$ ]]; then
            cluster_number="${BASH_REMATCH[1]}"
            group_number='unique'
            cluster_group="${cluster_number}_unique"
        elif [[ "$stem" =~ ^asinput_([^_]+)_([^_]+)$ ]]; then
            cluster_number="${BASH_REMATCH[1]}"
            group_number="${BASH_REMATCH[2]}"
            cluster_group="${cluster_number}_${group_number}"
        fi

        stats="$(awk '
            /^>/ {
                if (seen) {
                    n += 1
                    sum += sequence_length
                    sumsq += sequence_length * sequence_length
                }
                sequence_length = 0
                seen = 1
                next
            }
            {
                gsub(/[[:space:]]/, "", $0)
                sequence_length += length($0)
            }
            END {
                if (seen) {
                    n += 1
                    sum += sequence_length
                    sumsq += sequence_length * sequence_length
                }
                if (n > 0) {
                    average = sum / n
                    if (n > 1) {
                        variance = (sumsq - (sum * sum / n)) / (n - 1)
                        if (variance < 0 && variance > -0.0000001) {
                            variance = 0
                        }
                        standard_deviation = sqrt(variance)
                    } else {
                        standard_deviation = 0
                    }
                    printf "%d\t%.6f\t%.6f", n, average, standard_deviation
                } else {
                    printf "0\tNA\tNA"
                }
            }
        ' "$fasta_file")"

        kraken_record_type="${kraken_record_type_by_file[$asinput_file]:-none}"
        kraken_record_count="${kraken_record_count_by_file[$asinput_file]:-0}"
        kraken_classified_count="${kraken_classified_count_by_file[$asinput_file]:-0}"
        kraken_unclassified_count="${kraken_unclassified_count_by_file[$asinput_file]:-0}"
        kraken_taxids="${kraken_taxids_by_file[$asinput_file]:-none}"
        kraken_taxon_names="${kraken_taxon_names_by_file[$asinput_file]:-none}"
        kraken_taxid_counts="${kraken_taxid_counts_by_file[$asinput_file]:-none}"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$barcode" \
            "$asinput_file" \
            "$cluster_group" \
            "$cluster_number" \
            "$group_number" \
            "$stats" \
            "$kraken_record_type" \
            "$kraken_record_count" \
            "$kraken_classified_count" \
            "$kraken_unclassified_count" \
            "$kraken_taxids" \
            "$kraken_taxon_names" \
            "$kraken_taxid_counts" \
            "$fasta_file" >> "$OUTPUT_FILE"
        done
done

echo "Wrote $OUTPUT_FILE"
