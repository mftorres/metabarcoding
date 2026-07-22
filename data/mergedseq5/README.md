# `mergedseq5` barcode output tables

This directory contains per-barcode outputs from primer tracking, amplicon-sorter read/group tracing, Kraken taxonomic assignment, and the `evaluate_tagtracking_taxid_reduced.ipynb` notebook.

Most files are named with the barcode prefix:

```text
barcodeXX_<table_name>.tsv
```

For example:

```text
barcode01_top10_taxa_per_forward_tag.tsv
barcode02_target_species_detection.tsv
```

All tables documented here are tab-delimited, even when an upstream file may have used a `.csv` extension elsewhere in the workflow.

## Amplicon-sorter length-check table

### `mergedseq5_amplicon_sorter_length_check.tsv`

Mergedseq5-wide table produced by:

```text
code/length_check_amplicon_sorter.sh
```

Run the script after amplicon_sorter and Kraken2 have completed for the barcode folders. The script uses the per-barcode `filtered/combined/asinput_*_*.fasta` files plus the Kraken2 read-map and classification files:

```text
data/mergedseq5/barcodeXX/filtered/combined/asinput_consensussequences_noNUMT_200filter.read_map.tsv
data/mergedseq5/barcodeXX/filtered/combined/asinput_consensussequences_noNUMT_200filter.fasta.kraken
data/mergedseq5/barcodeXX/filtered/combined/asinput_consensussequences_noNUMT_200filter.fasta.report
```

The full `code/run_metabarcoding_pipeline.sh` step includes amplicon_sorter, Kraken2, and the Kraken read map, so the length check can be run after that pipeline finishes.

Default command:

```bash
code/length_check_amplicon_sorter.sh
```

Default output:

```text
data/mergedseq5/mergedseq5_amplicon_sorter_length_check.tsv
```

The script scans:

```text
data/mergedseq5/barcode*/filtered/combined/asinput_*_*.fasta
```

Kraken2 taxid assignment is attached as follows:

- For grouped files such as `asinput_1_0.fasta`, the script uses the corresponding consensus sequence in `asinput_consensussequences_noNUMT_200filter.fasta.kraken`. The link is made through `asinput_consensussequences_noNUMT_200filter.read_map.tsv`, where consensus read IDs follow the `consensus_combinedrawreads_<cluster>_<group>(...)` pattern.
- Unique-read FASTA files are not sent to Kraken2 in the current consensus-only pipeline, so they do not receive direct Kraken taxids in this table.
- Taxon names are taken from `asinput_consensussequences_noNUMT_200filter.fasta.report`.

Columns:

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

## Consensus noNUMT FASTA files

### `barcodeXX/filtered/combined/asinput_consensussequences_noNUMT_200filter.fasta`

Per-barcode FASTA file produced by the noNUMT step inside:

```text
code/run_metabarcoding_pipeline.sh
```

This step runs after amplicon_sorter has produced:

```text
barcodeXX/filtered/combined/asinput_consensussequences.fasta
```

The pipeline translates consensus sequences in all six frames with `seqkit translate --frame 6 --append-frame --trim`, then keeps amino-acid sequences at or above `NONUMT_MIN_AA_LENGTH` using `seqkit seq --min-len`.

The amino-acid output is used as the noNUMT detector. Its headers can contain `_frame=#` or `_frame=-#`; the pipeline strips that frame text, deduplicates passing consensus IDs, and extracts the corresponding original nucleotide consensus sequences. The cleaned nucleotide output is the file used as Kraken2 input.

The outputs are written beside the input consensus FASTA:

```text
barcodeXX/filtered/combined/asinput_consensussequences_aa_noNUMT_200filter.fasta
barcodeXX/filtered/combined/asinput_consensussequences_noNUMT_200filter.fasta
```

The configurable variables are in `config/pipeline.env`:

```text
NONUMT_MIN_AA_LENGTH
NONUMT_TRANSLATION_TABLE
NONUMT_SEQKIT_THREADS
SEQKIT_CONDA_ENV
```

`SEQKIT_CONDA_ENV` defaults to `umitracker`. If `seqkit` is not available in the active shell, the pipeline runs it with `conda run -n "$SEQKIT_CONDA_ENV" seqkit`.

## Amplicon-sorter primer-diversity table

### `mergedseq5_amplicon_sorter_primer_diversity.tsv`

Mergedseq5-wide table produced by:

```text
code/primer_diversity_amplicon_sorter.sh
```

Run the script after amplicon_sorter has produced the per-barcode `filtered/combined/asinput_*_*.fasta` files and after seqkit primer/UMI hits have been parsed into the combined best-hit table:

```text
data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_YYYYMMDD.txt
```

The script uses `filtered/combined/amplicon_sorter_read_map.tsv` to translate amplicon_sorter FASTA record IDs back to original ONT read IDs, then joins those ONT read IDs to the seqkit best-hit table.

Default command:

```bash
code/primer_diversity_amplicon_sorter.sh
```

The default command automatically uses the latest matching best-hit table under `data/primertrack/`. To force a specific dated best-hit table:

```bash
code/primer_diversity_amplicon_sorter.sh \
  data/mergedseq5 \
  data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_YYYYMMDD.txt
```

Default output:

```text
data/mergedseq5/mergedseq5_amplicon_sorter_primer_diversity.tsv
```

The script scans:

```text
data/mergedseq5/barcode*/filtered/combined/asinput_*_*.fasta
```

Columns:

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

## Read-level tracking tables

### `barcodeXX_read_primer_group_taxid_reduced.tsv`

Reduced read-level table for all raw ONT reads in a barcode. This is the main input used by the evaluation notebook as `tracking`.

| Column | Meaning |
|---|---|
| `raw_ont_read_id` | Original raw ONT read identifier. Read-count summaries deduplicate this column. |
| `forward_patternName` | Best detected forward-primer hit for the read, after primer-name conversion to PCR primer names. Empty means no forward primer was detected. |
| `reverse_patternName` | Best detected reverse-primer hit for the read, after primer-name conversion to PCR primer names. Empty means no reverse primer was detected. |
| `group_consensus_kraken_taxid` | Kraken taxid assigned to the amplicon-sorter group consensus. This value is propagated to reads that belong to that group. |

### `barcodeXX_insecta_read_primer_group_taxid_reduced.tsv`

Reduced read-level table for the Insecta subset only. This is the main input used by the evaluation notebook as `insecta_tracking`.

Columns are the same as `barcodeXX_read_primer_group_taxid_reduced.tsv`.

The table is generated by intersecting the full read tracking table with reads/groups represented in the parsed Insecta Kraken output, tracing amplicon-sorter identifiers back to original ONT read IDs.

The parsed Insecta Kraken output is a dated file, for example `data/kraken_parsed/mergedseq5_insecta_kraken_reads_YYYYMMDD.csv`. When a new `all_kraken_reads` or `insecta_kraken_reads` table is generated by `code/kraken_parser.ipynb`, update downstream code and commands to point to the new dated filename. In practice, this mainly means updating the `--insecta-kraken-reads` argument used by `code/build_barcode_read_tracking_table.py`. If an older dated file is reused, these Insecta reduced tables can be stale, empty for later barcodes, or inconsistent with the current Kraken parsing run.

### `barcodeXX_primer_hits.tsv`

Primer-hit table produced before best-hit reduction. This table records primer matches detected in barcode reads.

| Column | Meaning |
|---|---|
| `seqID` | Read or sequence identifier used in the primer-hit search. |
| `patternName` | Primer pattern name detected in the read. |
| `pattern` | Primer sequence or pattern searched. |
| `strand` | Strand/orientation of the detected hit. |
| `start` | Start coordinate of the primer hit in the read/sequence. |
| `end` | End coordinate of the primer hit in the read/sequence. |
| `matched` | Matched sequence segment. |

## Notebook summary/export tables

These tables are exported by `code/evaluate_tagtracking_taxid_reduced.ipynb`.

### `barcodeXX_primer_detection_classification_summary.tsv`

Summary of forward-primer detection and Kraken classification status for the complete tracking table and the Insecta-only table.

| Column | Meaning |
|---|---|
| `forward_primer_detection` | `detected` if at least one row for the raw read has a forward-primer hit; `not_detected` otherwise. |
| `classification` | `classified` if the read has a non-zero Kraken taxid; `unclassified` if taxid is missing or `0`. |
| `tracking_read_count` | Number of distinct raw ONT reads in the full tracking table for that primer/classification category. |
| `tracking_percentage` | Percentage within the corresponding `tracking` forward-primer-detection category. |
| `insecta_tracking_read_count` | Number of distinct raw ONT reads in the Insecta-only table for that primer/classification category. |
| `insecta_tracking_percentage` | Percentage within the corresponding `insecta_tracking` forward-primer-detection category. |

### `barcodeXX_top10_taxa_per_forward_tag.tsv`

One row per forward primer. Lists expected sample/species metadata and the ten most common classified Insecta taxids associated with that forward primer.

| Column | Meaning |
|---|---|
| `forward_patternName` | Forward primer name after conversion to PCR primer names. |
| `expected_forward` | Forward primer expected from the PCR design metadata for this sample/primer row. |
| `expected_reverse` | Reverse primer expected from the PCR design metadata for this sample/primer row. |
| `reverse_primer` | Backward-compatible copy of the expected reverse primer associated with the forward primer/sample metadata. |
| `sample_id` | SaproxDNA sample code associated with the primer/sample metadata. |
| `species` | Expected target species for the sample code. Multiple species are separated with `;` when a sample code maps to more than one species. |
| `species_genus` | Genus extracted from the expected target species name. Multiple values are separated with `;` when applicable. |
| `species_family` | Expected family resolved through NCBI taxonomy, using species first and genus fallback when the species name is absent from NCBI. |
| `target_genus_detected` | Whether any taxid associated with the forward primer matches the expected target genus. |
| `target_detected_reverse` | Reverse primer associated with the detected target-genus taxid for this forward primer. Empty if no target genus was detected. |
| `target_detected_readcount` | Distinct raw-read count supporting the detected target-genus taxid for this forward primer. Empty if no target genus was detected. |
| `taxon_1` | Scientific name of the most common classified Insecta taxid for the forward primer. |
| `taxon_1_genus` | Genus of `taxon_1`, resolved through NCBI taxonomy. |
| `taxon_1_family` | Family of `taxon_1`, resolved through NCBI taxonomy. |
| `read_count_1` | Distinct raw-read count supporting `taxon_1`. |
| `taxon_2` | Scientific name of the 2nd most common classified Insecta taxid for the forward primer. |
| `read_count_2` | Distinct raw-read count supporting `taxon_2`. |
| `taxon_3` | Scientific name of the 3rd most common classified Insecta taxid for the forward primer. |
| `read_count_3` | Distinct raw-read count supporting `taxon_3`. |
| `taxon_4` | Scientific name of the 4th most common classified Insecta taxid for the forward primer. |
| `read_count_4` | Distinct raw-read count supporting `taxon_4`. |
| `taxon_5` | Scientific name of the 5th most common classified Insecta taxid for the forward primer. |
| `read_count_5` | Distinct raw-read count supporting `taxon_5`. |
| `taxon_6` | Scientific name of the 6th most common classified Insecta taxid for the forward primer. |
| `read_count_6` | Distinct raw-read count supporting `taxon_6`. |
| `taxon_7` | Scientific name of the 7th most common classified Insecta taxid for the forward primer. |
| `read_count_7` | Distinct raw-read count supporting `taxon_7`. |
| `taxon_8` | Scientific name of the 8th most common classified Insecta taxid for the forward primer. |
| `read_count_8` | Distinct raw-read count supporting `taxon_8`. |
| `taxon_9` | Scientific name of the 9th most common classified Insecta taxid for the forward primer. |
| `read_count_9` | Distinct raw-read count supporting `taxon_9`. |
| `taxon_10` | Scientific name of the 10th most common classified Insecta taxid for the forward primer. |
| `read_count_10` | Distinct raw-read count supporting `taxon_10`. |

### `barcodeXX_unsorted_taxa.tsv`

Taxid summary for classified Insecta reads without a detected forward primer.

| Column | Meaning |
|---|---|
| `taxid` | Kraken/NCBI taxid detected among reads with no forward-primer hit. |
| `taxon_name` | Scientific name for `taxid`. |
| `taxid_genus` | Genus resolved for `taxid` through NCBI taxonomy. |
| `read_count` | Number of distinct raw ONT reads without a detected forward primer supporting this taxid. |
| `target_species` | Expected target species whose genus matches `taxid_genus`. Multiple matches are separated with `;`. Empty if no expected target genus matches. |
| `target_sxdcodes` | SaproxDNA sample codes corresponding to `target_species`. Order matches the order of `target_species`. |
| `target_forward_primers` | Forward primers associated with the matching target species. Multiple values are separated with `;`. |
| `target_reverse_primers` | Reverse primers associated with the matching target species. Multiple values are separated with `;`. |

### `barcodeXX_target_species_detection.tsv`

Target-focused table. Rows contain expected target species and the taxids whose genus matches each target genus, with and without detected forward primers.

| Column | Meaning |
|---|---|
| `target_species` | Expected target species from the sample/species metadata. |
| `target_genus` | Genus extracted from `target_species`. |
| `target_family` | Expected target family resolved through NCBI taxonomy, using species first and genus fallback. |
| `target_sxdcodes` | SaproxDNA sample code(s) for the target species. Multiple values are separated with `;`. |
| `expected_forward` | Forward primer(s) expected from the PCR design metadata in `PCR_to_species_*`, using the selected PCR code/sample code. Multiple values are separated with `;`. |
| `expected_reverse` | Reverse primer(s) expected from the PCR design metadata in `PCR_to_species_*`, using the selected PCR code/sample code. Multiple values are separated with `;`. |
| `forward_primer_detected` | `yes` if the supporting read had a detected forward primer; `no` if the supporting read had no detected forward primer. |
| `forward_primer` | Seqkit-detected forward primer for the row. Empty when `forward_primer_detected` is `no` or when the target species was not detected. |
| `forward_primer_agrees_expected` | `yes` if the seqkit-detected forward primer is present in `expected_forward`; `no` if it is not. Empty when no forward primer was detected. |
| `detected_reverse_primers` | Seqkit-detected reverse primer(s) among reads supporting the row. Multiple values are separated with `;`. |
| `reverse_primer_agrees_expected` | `yes` if at least one seqkit-detected reverse primer is present in `expected_reverse`; `no` if none are. Empty when no reverse primer was detected. |
| `taxid_detected` | Kraken/NCBI taxid whose genus matches `target_genus`. Empty when the expected target species was not detected. |
| `taxid_primer_read_count` | Distinct raw-read count for this exact target species, forward-primer status/primer, and detected taxid combination. |
| `taxid_read_count` | Distinct raw-read count for this target species and detected taxid across primer statuses. |
| `forward_primer_read_count` | Distinct raw-read count for this target species and detected forward primer. Empty for no-forward-primer rows. |
| `no_forward_primer_read_count` | Distinct raw-read count for this target species among reads with no detected forward primer. |
| `target_species_read_count` | Total distinct raw-read count supporting the target species by target-genus match. |
| `taxid_detected_name` | Scientific name of `taxid_detected`. |
| `taxid_detected_genus` | Genus resolved for `taxid_detected` through NCBI taxonomy. |
| `taxid_in_50pct_forward_primers` | `yes` if `taxid_detected` is a widespread taxid found in at least 50% of observed forward primers; `no` otherwise. Empty when no taxid was detected. |

### `barcodeXX_taxids_across_50pct_forward_primers.tsv`

Taxids detected with at least 50% of observed forward primers in the Insecta tracking table. These are treated as widespread taxa because they occur across many primer/sample combinations.

| Column | Meaning |
|---|---|
| `taxid` | Kraken/NCBI taxid. |
| `taxid_name` | Scientific name for `taxid`. |
| `taxid_genus` | Genus resolved for `taxid` through NCBI taxonomy. |
| `forward_primer_count` | Number of distinct forward primers where this taxid was detected. |
| `forward_primer_percentage` | Percentage of observed forward primers where this taxid was detected. |
| `forward_primer_detected` | Semicolon-separated list of forward primers where this taxid was detected. |
| `forward_primer_not_detected` | Semicolon-separated list of observed forward primers where this taxid was not detected. Empty if detected in all observed forward primers. |
| `expected_forward` | Expected forward primer(s), from PCR design metadata, corresponding to the detected forward primers for this taxid. |
| `expected_reverse` | Expected reverse primer(s), from PCR design metadata, corresponding to the detected forward primers for this taxid. |
| `forward_primer_agrees_expected` | `yes` if at least one detected forward primer is present in `expected_forward`; `no` otherwise. |
| `detected_reverse_primers` | Seqkit-detected reverse primer(s) among reads supporting this taxid. Multiple values are separated with `;`. |
| `reverse_primer_agrees_expected` | `yes` if at least one detected reverse primer is present in `expected_reverse`; `no` otherwise. |
| `forward_primer_read_count` | Total distinct raw-read count supporting this taxid among reads with detected forward primers. |
| `detected_without_forward_primer` | `yes` if the same taxid is also detected among reads with no forward-primer hit; `no` otherwise. |
| `no_forward_primer_read_count` | Distinct raw-read count supporting this taxid among reads with no detected forward primer. |

### `barcodeXX_target_species_taxids_below_50pct_forward_primers.tsv`

Filtered version of `barcodeXX_target_species_detection.tsv`. It keeps target-species detections whose detected taxid is not widespread across at least 50% of forward primers.

Columns are the same as `barcodeXX_target_species_detection.tsv`.

Interpretation: this table is useful for focusing on target-genus detections that are less likely to be broad, cross-primer background taxa.

### `barcodeXX_target_species_detection_trusted.tsv`

Expanded and trust-annotated version of `barcodeXX_target_species_detection.tsv`, produced by `code/sankey_target_species_detection.ipynb`.

The table keeps the original input columns and adds one row per expanded expected-forward-primer/detected-forward-primer combination. It is the input used for both Sankey figures.

Important added columns:

| Column | Meaning |
|---|---|
| `source_row_index` | Zero-based row index from the original `barcodeXX_target_species_detection.tsv` table. |
| `target_species_original` | Original `target_species` cell before expansion. |
| `target_sxdcodes_original` | Original semicolon-separated `target_sxdcodes` cell before expansion. |
| `expected_forward_original` | Original semicolon-separated `expected_forward` cell before expansion. |
| `forward_primer_original` | Original detected `forward_primer` cell before expansion. |
| `target_sxdcode` | Single expanded SaproxDNA code. When `target_sxdcodes` and `expected_forward` lists are the same length, this remains positionally paired with `expected_forward_primer`. |
| `expected_forward_primer` | Single expanded expected forward primer. |
| `detected_forward_primer` | Single expanded detected forward primer from seqkit. Missing detected primers are `none`. |
| `detected_species` | Detected species/scientific name from `taxid_detected_name`. Missing detected species are `none`. |
| `detected_genus` | Detected genus from `taxid_detected_genus`. |
| `target_genus_clean` | Cleaned target genus used for trust comparison. |
| `expansion_note` | How `target_sxdcodes` and `expected_forward` were expanded: positional pairing, repeated single value, or all combinations after length mismatch. |
| `sankey_read_count` | Read-count weight used in the Sankey after splitting the original row's count across expanded combinations. |
| `target_expected_pair_key` | Internal normalized key for the expanded `target_species + expected_forward_primer` pair. |
| `detected_species_key` | Internal normalized detected-species key used for uniqueness checks. |
| `target_genus_key` | Internal normalized target-genus key. |
| `detected_genus_key` | Internal normalized detected-genus key. |
| `trusted` | Boolean trust classification. |
| `trust_reason` | Rule that determined the trust classification. |

Trust rules:

- Primer-detected rows are trusted only when `expected_forward_primer == detected_forward_primer`.
- For a `target_species + expected_forward_primer` pair with multiple detected forward primers, only rows where the detected primer matches the expected primer are trusted.
- For target species with multiple expected forward primers, each expanded expected/detected primer row is evaluated independently.
- Rows with no detected forward primer are trusted only when `target_genus` and `taxid_detected_genus` match and the detected species is unique to one `target_species + expected_forward_primer` pair in that barcode.
- Rows with no detected forward primer are not trusted when the detected species appears in multiple target/species expected-primer pairs.

### `mergedseq5_target_species_detection_trusted.tsv`

Combined trusted table exported by the final section of `code/sankey_target_species_detection.ipynb`.

Path:

```text
data/mergedseq5/mergedseq5_target_species_detection_trusted.tsv
```

This table concatenates trusted rows from all selected barcode-level `barcodeXX_target_species_detection_trusted.tsv` tables. If a barcode-level trusted table does not exist when the notebook is run, the notebook rebuilds it in memory from `barcodeXX_target_species_detection.tsv` using the same expansion and trust rules.

Additional columns relative to the barcode-level trusted table:

| Column | Meaning |
|---|---|
| `barcode` | Barcode folder from which the trusted row came. |
| `source_target_species_detection_file` | Source `barcodeXX_target_species_detection.tsv` file used for the row. |

### `barcodeXX_expected_species_not_detected.tsv`

Expected target species for which no matching target-genus taxid was detected in the Insecta tracking table.

| Column | Meaning |
|---|---|
| `target_species` | Expected target species from the sample/species metadata. |
| `target_genus` | Genus extracted from `target_species`. |
| `target_family` | Expected target family resolved through NCBI taxonomy, using species first and genus fallback. |
| `target_sxdcodes` | SaproxDNA sample code(s) for the expected species. Multiple values are separated with `;`. |
| `expected_forward` | Forward primer(s) expected for the target species/sample metadata. Multiple values are separated with `;`. |
| `expected_reverse` | Reverse primer(s) expected for the target species/sample metadata. Multiple values are separated with `;`. |

## Figure exports

The evaluation notebook also writes PNG figures at 250 dpi under:

```text
barcodeXX/figures/
```

Current figure names are:

```text
barcodeXX_primer_detection_classification.png
barcodeXX_reads_per_forward_primer.png
barcodeXX_reads_per_reverse_primer.png
barcodeXX_unique_taxids_per_forward_primer.png
barcodeXX_insecta_reads_vs_taxid_diversity.png
barcodeXX_top10_taxids_per_forward_primer.png
```

Additional Sankey/alluvial figures can be generated with:

```text
code/sankey_target_species_detection.ipynb
```

This notebook reads the full `barcodeXX_target_species_detection.tsv` table for each barcode and writes an expanded trust-annotated table:

```text
barcodeXX_target_species_detection_trusted.tsv
```

It also writes two Sankey figures:

```text
barcodeXX/figures/barcodeXX_target_detection_sankey.png
barcodeXX/figures/barcodeXX_target_detection_sankey_trusted.png
```

The same notebook also reads the filtered below-50% table:

```text
barcodeXX_target_species_taxids_below_50pct_forward_primers.tsv
```

For that filtered table, it does not export a new expanded/trusted TSV. It only writes Sankey figures:

```text
barcodeXX/figures/barcodeXX_target_species_taxids_below_50pct_forward_primers_sankey.png
barcodeXX/figures/barcodeXX_target_species_taxids_below_50pct_forward_primers_sankey_trusted.png
```

Sankey rules:

- The diagram flows from `expected_forward_primer` to `detected_forward_primer` to `detected_species` to `target_species`.
- `expected_forward_primer` is taken directly from the row's `expected_forward` column; missing or empty values are grouped as `none`. Multiple semicolon-separated expected primers are split.
- `target_species` is taken directly from the same row's `target_species` column; it is not globally re-derived from the expected primer. Missing or empty values are grouped as `none`.
- `detected_species` is taken from `taxid_detected_name`; missing or empty values are grouped as `none`.
- `detected_forward_primer` is taken from `forward_primer`; missing or empty detected forward primers are grouped as `none`.
- `target_sxdcodes` and `expected_forward` are expanded together. If both contain the same number of semicolon-separated values, they are paired by position so the first SxD code remains linked to the first expected forward primer, and so on. If one side has one value and the other has multiple, the single value is repeated. Otherwise, all combinations are used and the row is flagged in `expansion_note`.
- Link weights are based on read-count columns from `barcodeXX_target_species_detection.tsv`: `taxid_primer_read_count` for detected-primer rows, `no_forward_primer_read_count` for no-forward-primer rows, with fallback to `taxid_read_count` and then `target_species_read_count`.
- If one row expands to multiple expected-primer/detected-primer/species combinations, the row's read count is split evenly across those combinations as `sankey_read_count` so the total flow is conserved.
- Detected species and target species blocks are sorted alphabetically, with `none` at the end.
- Expected and detected forward-primer blocks are initialized alphanumerically/naturally, then refined with weighted barycentric sweeps to reduce path overlap/crossings. This is a layout-only operation and does not change link weights.
- The Sankey plot is exploratory. It shows expected-primer, expected-target, detected-taxid, and detected-primer associations, but it does not by itself prove that the detected taxid is the exact expected species.

Trust rules used for `barcodeXX_target_species_detection_trusted.tsv` and the green trusted Sankey:

- After expansion, primer-detected rows are trusted only when `expected_forward_primer` and `detected_forward_primer` are the same.
- If a `target_species + expected_forward_primer` pair has multiple detected forward primers, only rows where the detected primer matches the expected primer are trusted; mismatching detected primers are not trusted.
- If a target species has multiple expected forward primers, each expanded expected/detected primer row is evaluated independently.
- Rows with no detected forward primer are trusted only when `target_genus` and `taxid_detected_genus` match and the detected species is unique to one `target_species + expected_forward_primer` pair in that barcode.
- Rows with no detected forward primer are not trusted when the detected species appears in multiple target/species expected-primer pairs.
- The trusted Sankey colors trusted links in green and untrusted links in grey.

## Notes on counts and classification

- Read counts are distinct counts of `raw_ont_read_id`, not row counts.
- Kraken taxid `0` is treated as unclassified.
- In grouped amplicon-sorter reads, `group_consensus_kraken_taxid` is propagated to reads in the group.
- Unique-read FASTA files are not sent to Kraken2 in the current consensus-only pipeline.
- The evaluation notebook focuses most taxonomy summaries on `insecta_tracking`, while some primer-detection comparisons include both the complete `tracking` table and `insecta_tracking`.
