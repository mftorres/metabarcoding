# Module 2: Phylogenetics

[Module 1](pipeline_from_fastplong_onwards.md) cleans reads, assigns taxonomy with
SINTAX and generates read-level analyses. Module 2 exports selected cleaned reads,
obtains external references, builds cluster consensus sequences, optionally screens
translations, aligns nucleotides, estimates an IQ-TREE phylogeny and draws a Baltic
midpoint-rooted divergence tree. Run from the repository root in `taxonomy`.

## Run the module

```bash
conda activate taxonomy
python code/postprocess_sintax.py --bold \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_bold_run02 --threads 8
```

Use a new output directory for a new analysis. Dependencies are `vsearch`, `seqkit`,
`iqtree`, Python with Baltic >=1.0.0 and Matplotlib, and the chosen aligner: `mafft`
(default) or `macse`. NCBI mode additionally requires Entrez Direct `esearch` and
`efetch`. MACSE is available in the `taxonomy` environment (tested with 2.07).

## Choose the reference dataset

Exactly one source flag is required:

- `--download`: resolve a name or numeric taxid with NCBI, then download COI/CO1/COX1
  nucleotide records including descendants. Searches also include cytochrome
  oxidase I/1 and cytochrome c oxidase subunit I/1 title aliases.
- `--bold`: stream `.fasta.gz`, `.fa.gz` and `.fna.gz` databases in `data/bold/`
  locally. `--bold-db PATH` selects one file or another directory. This mode accepts
  a taxon name, not a numeric taxid, and makes no Entrez calls. It matches a substring
  of the full FASTA header, ignoring case and treating underscores as spaces.
  Every matching marker is included; it does not implicitly restrict markers to COI.

BOLD headers are expected to contain `processID|marker|country|comma-separated lineage`.
The original header stays in the FASTA description. Safe IDs `BOLD_000000001`, etc.,
are mapped to process IDs and taxonomy by `bold_reference_metadata.tsv` with columns
`sequence_id`, `process_id`, `marker`, `taxon`, `original_header`, `source_file`.
The display taxon is the deepest nonempty lineage entry.

Mitogenomes are excluded by title by default; `--include-mitogenomes` permits them.
This title-based rule is not a guarantee that every mitogenome is detected.
Both sources are then filtered with `seqkit seq --min-len 500 --max-len 2000`:
the bounds are inclusive. This reference length filter remains active even when
translation is skipped. Empty reference datasets stop the run.

## Local read selection and headers

`export_sintax_taxa.py` selects reads when the queried name, complete MIDORI label
or taxid matches a rank in the full SINTAX lineage **and that rank has confidence
>=0.80**. Lower ranks need not pass. Tipulidae at 0.90 with Tipula at 0.20 qualifies
for a Tipulidae query but not a Tipula query. Missing confidence does not qualify.
SINTAX confidence is bootstrap support, not a calibrated probability.

Cleaned sequences are extracted by read ID. Raw reads are not used as phylogenetic
input. Headers contain:

```text
folder-forwardprimer-reverseprimer-lowest_trusted_taxon-contamination_flag-readID
```

The taxon is the lowest rank with confidence >=0.80; contamination applies to that
trusted taxon according to the completed analysis. Primer names use the canonical
mapping in that analysis. Header components escape delimiter characters; read IDs
retain their tracking information. Extraction counts for every requested folder
and the total are saved in `extraction_counts.tsv` and `extraction.log`.

## Clustering: 98% consensus sequences

Local reads are partitioned by the exact first **four** header fields:
**(folder/barcode, forward primer, reverse primer, trusted taxonomic label)**.
All reads in one tuple are submitted together to one VSEARCH invocation. Different
tuples never share a cluster. Contamination flag and read ID do not create extra
groups. `none` and `unknown` remain distinct primer values.

External references are clustered globally, separately from local reads. Keeping
these datasets separate prevents local reads and reference records from being
mixed in the same cluster. Both use:

```bash
vsearch --cluster_fast input.fasta --id 0.98 --iddef 2 \
  --strand both --qmask none --notrunclabels \
  --consout consensus.fasta --uc clusters.uc --threads 8
```

The module sets minimum length 1 and maximum length to the longest input so this
step adds no length filter. Identity definition 2 excludes terminal gaps.
Downstream sequences are **consensuses**, potentially different from individual
input reads. VSEARCH's `centroid=` and `;seqs=N` label decorations are removed;
the seed header is retained for tracking. Its taxon, contamination and read ID or
accession describe the **seed**, not a new classification of the consensus.
No new SINTAX classification is performed. UC files identify all cluster members.

Outputs include:

- `cleaned_TAXON_DATE.vsearch98.fasta`: combined local consensus sequences.
- `COI_TAXON_DATE.vsearch98.fasta` or `BOLD_TAXON_DATE.vsearch98.fasta`: reference consensuses.
- `vsearch_local_groups.tsv`: `group_id`, `source_folder`, `forward_primer`,
  `reverse_primer`, `taxon_label`, `input_sequences`, `consensus_sequences`, `identity`.
- `vsearch_local_groups/`: each group's inputs, consensus FASTA, UC and stdout log.
- `vsearch_local_memberships.tsv`: group ID followed by the ten UC fields.
  Cluster numbers are relative to their group.
- Reference `.vsearch98.uc` and `vsearch_references.stdout.log`.
- `clustering_counts.tsv`: `dataset`, `input_sequences`, `consensus_sequences`, `identity`.

The local and reference consensus sets are concatenated to `combined_TAXON_DATE.fasta`.
There is no further clustering after combining them.

## Translation selection and nuclear sequences

By default, SeqKit translates the combined consensuses in all six frames with
`--translation-table 5` (invertebrate mitochondrial). At least one frame must
contain **at least 200 consecutive amino acids**, splitting at stops `*` and
unknown residues `X`. `--min-aa` adjusts the minimum. Outputs are
`combined_TAXON_DATE.six_frames.faa`, `translation_lengths.tsv`, and the retained
nucleotides `combined_TAXON_nonumt_DATE.fasta`. This is a screening heuristic,
not proof of the absence of NUMTs. The original consensus nucleotides are retained.

**`--skip-translation`** bypasses both translation and amino-acid selection and
keeps all consensuses in `combined_TAXON_selected_DATE.fasta`. It creates no
six-frame FASTA or translation-length table. Use this with default MAFFT for
**nuclear loci, particularly noncoding sequences**, for which the mitochondrial
translation screen is inappropriate:

```bash
python code/postprocess_sintax.py --bold --bold-db data/bold/nuclear_marker.fasta.gz \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_nuclear_run01 \
  --skip-translation --threads 8
```

The nuclear database path above is a placeholder for your marker-specific database,
using the BOLD header format. Skipping translation does not change marker selection:
**NCBI `--download` still searches COI**, and BOLD mode includes all matching markers.
The local reads, their SINTAX reference assignments, and the external database must
be appropriate to the locus being analysed. The 500–2,000 nt reference filter remains.

## Alignment choice

MAFFT is the default (`--mafft` may also be explicit). Current arguments are:

```bash
mafft --retree 2 --maxiterate 0 --adjustdirectionaccurately \
  --op 20 --ep 18 --thread 8 selected_consensus.fasta
```

This is FFT-NS-2 nucleotide alignment, suitable for nuclear sequences as well as
other nucleotide loci. Strand-adjusted headers are restored and recorded in
`alignment_orientation.tsv`. All representative IDs, equal alignment lengths and
preservation of nucleotide sequences (allowing reverse complementation) are checked.

Add **`--macse`** for protein-coding alignment followed by end trimming:

```bash
python code/postprocess_sintax.py --bold --macse \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_macse_run01 --threads 8
```

The equivalent MACSE operations are:

```bash
macse -prog alignSequences -gc_def 5 -gap_ext 2.0 -gap_op 10 \
  -out_AA PREFIX.aa.aln.fasta -out_NT PREFIX.macse.aln.fasta \
  -seq selected_consensus.fasta
macse -prog trimAlignment -min_percent_NT_at_ends 0.7 \
  -out_NT PREFIX.macse.trimmed.aln.fasta \
  -out_trim_info PREFIX_trimmed.csv -align PREFIX.macse.aln.fasta
```

**MACSE 2.07 requires `0.7` for 70%; literal `70` fails.** `--translation-table`
sets the genetic code for both SeqKit and MACSE (default 5); nuclear coding sequences
normally need the appropriate nuclear code instead. MACSE is for coding sequences,
not general noncoding nuclear alignment. Input orientation must be appropriate for
coding alignment; this branch does not invoke MAFFT strand adjustment.

MACSE runs with temporary safe IDs, then restores all original headers. Its raw
AA, nucleotide and trimmed nucleotide alignments are retained. The trimmed
nucleotides feed the next steps as `PREFIX.aln.fasta`, with frameshift symbols `!`
replaced by `N` (missing/unknown nucleotide) for IQ-TREE compatibility. No columns
are removed by that replacement. `.frameshifts.tsv` records replacement counts;
`.macse_ids.tsv` maps temporary IDs to original headers. The native MACSE trim
report is CSV as requested; other pipeline tables are tab-delimited with headers.
See [MACSE alignment/export documentation](https://www.agap-ge2pop.org/exportalignment/).

`PREFIX` is `combined_TAXON_nonumt_DATE` when translation selection is active,
or `combined_TAXON_selected_DATE` with `--skip-translation`. MACSE uses the same
selected input as MAFFT. Genetic code 5 is not used when translation is skipped
and MAFFT is chosen.

## Phylogeny and figure

IQ-TREE consumes the final nucleotide `.aln.fasta`, with `-m MFP -b 1000`, seed 666
and the requested thread count. `--bootstrap-method ultrafast` uses `-B 1000`.
At least four retained consensus sequences are required. Outputs use `PREFIX`.

Baltic renders a midpoint-rooted divergence tree with branch lengths in
substitutions/site and its scale bar. Existing IQ-TREE support labels are displayed
on the corresponding splits, including slash-separated values. If `.treefile` has
no support and `.contree` does, the supported consensus topology is used. Supports
are preserved during rerooting; none are invented. Midpoint rooting does not
establish the biological root.

Tip labels are trusted taxon + read ID + folder for local cluster seeds; GenBank
organism + accession.version for downloaded seeds; or available BOLD taxon + process
ID for BOLD seeds. GenBank metadata is fetched once and cached as
`reference_records.gb.xml`; concatenated Entrez XML batches are parsed locally.
BOLD metadata requires no network.

Figure outputs are `PREFIX.figure.pdf`, `.svg`, `.midpoint.nwk`, `.tips.tsv` and
`.plot.json`. The tips table columns are `tree_tip_id`, `source`, `taxon_or_organism`,
`accession`, `readid`, `source_folder`, `display_label`; BOLD process IDs occupy
`accession` for BOLD rows. To plot independently use `code/plot_phylogeny.py --tree
PATH --analysis-dir PATH --reference-metadata PATH --output-prefix PATH`, with
GenBank XML or the BOLD metadata TSV as the metadata input.

## Resume and audit

Add `--resume` to the original command with the same output directory and settings.
Completed retained/selected sequences restart at alignment. A successfully completed
alignment restarts at IQ-TREE without repeating extraction or reference retrieval.
Completion hashes distinguish finished alignments from interrupted outputs. MAFFT
partial files are removed on interruption; MACSE works in temporary files and
publishes its completion marker only after both alignment and trimming validate.
Incomplete outputs are rebuilt on resume, with stale files archived under
`resume_archive/` and successful reference retrieval reused when needed.

Source, input database, taxon, folders, translation choice, aligner, genetic code
and other biological settings must match the saved run. Changed policies (including
98% four-field grouping) require a **new output directory**. Existing biological
outputs are not rewritten by a code update.

`parameters.json` stores settings, policy versions and executable paths;
`checkpoints.json` stores output hashes. `postprocessing.log` records commands,
stderr and counts. `postprocessing.failed` records failure; `postprocessing.complete`
is written only after phylogeny and plotting finish. DATE defaults to the run date
and may be specified with `--date YYYYMMDD`.

## Validation

Tests use temporary fixtures, never existing biological analyses:

```bash
conda activate taxonomy
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_postprocess_sintax.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_vsearch_groups.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_macse_alignment.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_bold_references.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_plot_phylogeny.py -v
```

## MACSE memory and alignment failures

The Conda MACSE launcher defaults to a 2 GB Java heap. Set `--macse-memory 32g`
(or another positive size ending in `m` or `g`) for larger alignments, leaving
sufficient RAM for the operating system and other processes. The module explicitly
passes this limit to both MACSE commands. The default remains `2g`.

Memory is a resource setting and may be increased during `--resume` in the same
output directory. Completed extraction, clustering and translation selection are
reused. For the Tipulidae BOLD run:

```bash
conda activate taxonomy
python code/postprocess_sintax.py --bold --macse --macse-memory 32g \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_bold_macse_run01 \
  --translation-table 5 --threads 8 --resume
```

Inspect `PREFIX.macse_align.stdout.log`, `PREFIX.macse_trim.stdout.log` and
`postprocessing.log` after a failure. A nonzero exit alone does not establish memory
exhaustion; inspect diagnostics before attributing a failure to memory. Resource
changes are saved in `resume_parameters.json`. Failed-attempt logs are archived
when retrying the alignment stage.
