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
`clipkit`, `iqtree`, Python with Baltic >=1.0.0 and Matplotlib, and the chosen aligner: `mafft`
(default) or `kc-align`. NCBI mode additionally requires Entrez Direct `esearch` and
`efetch`. kc-align is installed in the `taxonomy` environment.

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

`none` means a primer-tracking record exists for the read's raw parent, but no
primer was detected in that direction. `unknown` means no matching primer-tracking
record was found, so the primer assignment is unavailable. A missing record is
not evidence that the primer was absent. These meanings apply separately to the
forward and reverse primer fields.

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

**Important:** For now, the pipeline lacks a step for primer sequence removal. Before tree estimation, the user **must** remove the primer sequences from the alignment. In a next version of the pipeline, I will add the primer sequences as "adapter/barcode" to `fastplong` to remove them at early stages. I doubt (but will test) taxonomic classification will be severely affected by the presence of primer sequences in the clean reads (2026-Oct-29).

Add **`--kcalign --kcalign-reference PATH`** for reference-guided mixed-mode
codon alignment. The reference file is required: **provide a nucleotide sequence
you trust for this locus**, in the correct orientation and reading frame. It must
contain exactly one ungapped A/C/G/T sequence. Its length must be divisible by
three, and it must start at the first nucleotide of the first codon in the
reference coding region. Remove incomplete terminal codons using the known
reading frame; do not trim arbitrarily just to make the length divisible by three.
The reference must have no internal stops under the chosen genetic code. Automatic validation cannot establish
that the reference's taxonomic identity or biological annotation is correct.

**Important:** The reference must not contain COI primer sequences. kcalign looks for the reference pattern in each sequence and if the reference does not contain the primer sequences, primer sequences remaining in the clean reads are effectively removed before alignment and tree estimation.

```bash
conda activate taxonomy
python code/postprocess_sintax.py --bold --kcalign \
  --kcalign-reference data/postprocessing/coi_refs/SATIP962-10_Nigrotipula_nigra.fasta \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_bold_kcalign_run01 \
  --translation-table 5 --threads 6
```

The example requires a valid reference FASTA with a `>` header and a trusted,
correctly framed coding sequence. The module does not guess its frame or silently
trim bases. Divisibility by three and absence of internal stops alone do not prove
the biological reading frame; the user must verify its codon boundaries.

The alignment takes `combined_TAXON_nonumt_DATE.fasta`, the nucleotide consensus
sequences retained by the >=200 consecutive amino-acid screen. With
`--skip-translation`, it instead takes `combined_TAXON_selected_DATE.fasta`.
The corresponding prefix is `combined_TAXON_nonumt_DATE` or
`combined_TAXON_selected_DATE`. The mixed-mode command is equivalent to:

```bash
kc-align --mode mixed --reference trusted_reference.fasta \
  --sequences combined_TAXON_nonumt_DATE.fasta \
  --table 5 --keep --threads 6
```

**Installed-tool compatibility:** `--keep` is a switch, not a filename argument.
The tool writes `pre_align.fasta` and `codon_align.fasta` in its working directory.
Module 2 isolates those files in a temporary directory, restores tracking headers,
and publishes `PREFIX.pre_align.aa.fasta` and `PREFIX.aln.fasta`. The installed
CLI omits the threads argument when calling mixed mode; `code/kcalign_cli.py`
forwards it correctly to the installed library without changing the conda package.
The local parallel adapter also fixes an installed rejection-handling bug that
otherwise crashes when a worker rejects a read. It runs the original serial
per-read checks in isolated worker directories, preserves rejected IDs, and uses
a persistent pool with at most `threads / 3` workers. Biological checks remain unchanged.
Mixed mode requires **1 thread or a positive multiple of 3**. The module passes
its `--threads` value unchanged and rejects incompatible values such as 4 or 8.
Kalign and MAFFT (the library's fallback aligner) must also be available in `taxonomy`.

The trusted reference guides alignment but is not added as an extra phylogenetic
tip. An independent input consensus with the same sequence remains eligible.
Mixed mode can reject reads with suspected frameshifts, early stops or excessive
ambiguity even after the earlier translation screen. `PREFIX.kcalign_selection.tsv`
records every input ID as `aligned` or `rejected_by_kcalign`; details are in
`PREFIX.kcalign.stdout.log` and `postprocessing.log`.
`PREFIX.kcalign_ids.tsv` maps temporary tool IDs to original headers.
At least four aligned consensus sequences are required for IQ-TREE. ClipKIT trims
the published **nucleotide** `PREFIX.aln.fasta` before IQ-TREE; the retained pre-alignment AA file
is an audit output. ClipKIT smart-gap column trimming follows alignment as described below.

MAFFT remains the default for general nucleotide alignment. kc-align is for coding
sequences; skipping the earlier SeqKit filter does not disable kc-align's own
coding-sequence checks. `--translation-table` controls its genetic code (default 5).

## ClipKIT smart-gap trimming

After either MAFFT or kc-align, Module 2 automatically runs:

```bash
clipkit PREFIX.aln.fasta -m smart-gap -s nt -if fasta -of fasta \
  -o PREFIX.clipkit.aln.fasta
```

Smart-gap determines a gap threshold and removes gappy **alignment columns**.
It does not strip every gap independently from each sequence, and some gaps may
remain. The original alignment is preserved. Headers and the sequence set must
remain unchanged, and the trimmed alignment must be nonempty and rectangular.
IQ-TREE receives `PREFIX.clipkit.aln.fasta`. `clipkit.stdout.log` records the tool's
report; `clipkit.complete.json` records mode and input hash for resume validation.

This step is active for both aligners and with or without the translation screen.
No extra flag is required. On `--resume`, a completed untrimmed alignment can be
reused: missing or stale ClipKIT output is regenerated, and the dependent tree and
figure are rebuilt. Existing dependent outputs are archived. Failed partial trim
outputs are removed; extraction, clustering and alignment need not be repeated.

## Phylogeny and figure

For new runs, reference selection appends the source taxon label to the sequence
ID before length filtering, clustering and alignment: for example,
`BOLD_000000123__taxon__Tipula_maxima` or `AB123456.1__taxon__Tipula_maxima`.
Whitespace becomes underscores; punctuation unsafe for alignment/Newick tools
also becomes underscores. BOLD uses the deepest available taxon in its metadata
(which may be a genus or higher rank when species is unavailable); NCBI uses the
GenBank organism name. Original descriptions remain in the source files and audit
metadata, while selected reference headers contain only the annotated ID.
These names carry through both untrimmed and trimmed alignments into IQ-TREE.
Plotting accepts both original reference IDs and IDs with the `__taxon__` suffix,
using the original ID to look up source metadata. Local read headers are unchanged.
Resuming an older run reuses its existing sequence headers; this change does not
automatically rename completed runs.

IQ-TREE consumes the trimmed nucleotide `.clipkit.aln.fasta`, with `-m MFP`, seed 666
and the requested thread count. Bootstrap replicates are disabled; the pipeline
now estimates a maximum-likelihood tree without bootstrap branch support. The
1,000-replicate command remains commented in the source for reference. The
`--bootstrap-method` option is no longer accepted.
At least four retained consensus sequences are required. Outputs use `PREFIX`.

Baltic renders a midpoint-rooted divergence tree with branch lengths in
substitutions/site and its scale bar. Existing IQ-TREE support labels are displayed
on the corresponding splits, including slash-separated values. If `.treefile` has
no support and `.contree` does, the supported consensus topology is used. Supports
are preserved during rerooting; none are invented. Midpoint rooting does not
establish the biological root.

Tip labels are trusted taxon + read ID + folder + forward primer for local cluster seeds; GenBank
organism + accession.version for downloaded seeds; or available BOLD taxon + process
ID for BOLD seeds. GenBank metadata is fetched once and cached as
`reference_records.gb.xml`; concatenated Entrez XML batches are parsed locally.
BOLD metadata requires no network.

For the forward-primer field in local tip labels, `none` means tracking was
available but no forward primer was detected; `unknown` means the raw parent had
no matching tracking record. The label comes from the local cluster seed's
exported primer annotation. These values remain separate during clustering, as
described above; BOLD and GenBank tip labels do not include this field.

Figure outputs are `PREFIX.figure.pdf`, `.png`, `.svg`, `.midpoint.nwk`, `.tips.tsv` and
`.plot.json`. The tips table columns are `tree_tip_id`, `source`, `taxon_or_organism`,
`accession`, `readid`, `source_folder`, `display_label`; BOLD process IDs occupy
`accession` for BOLD rows. To plot independently use `code/plot_phylogeny.py --tree
PATH --analysis-dir PATH --reference-metadata PATH --output-prefix PATH`, with
GenBank XML or the BOLD metadata TSV as the metadata input.

## Resume and audit

Add `--resume` to the original command with the same output directory and settings.
Completed retained/selected sequences restart at alignment. A successfully completed
alignment restarts at ClipKIT, or at IQ-TREE if trimming is already complete,
without repeating extraction or reference retrieval.
Completion hashes distinguish finished alignments from interrupted outputs. MAFFT
partial files are removed on interruption; kc-align works in temporary files and
publishes its completion marker only after the nucleotide alignment validates.
Incomplete outputs are rebuilt on resume, with stale files archived under
`resume_archive/` and successful reference retrieval reused when needed.

Source, input database, taxon, folders, translation choice, aligner, genetic code
the trusted reference contents, and other biological settings must match the saved run. To change aligner or trusted reference while reusing validated retained sequences,
use `--resume --restart-alignment`. This archives the previous settings and rebuilds
alignment, tree and figure; subsequent resumes use the new settings. Other biological
settings must still match. Changed policies (including
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
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_kcalign_alignment.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_clipkit_alignment.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_bold_references.py -v
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_plot_phylogeny.py -v
```

On the first `--resume` after disabling bootstrap, completed alignment and ClipKIT
outputs are reused. The tree and figure are regenerated without bootstrap;
previous IQ-TREE outputs, including its checkpoint, are archived. A new
`phylogeny.complete.json` records the no-bootstrap setting and alignment hash.
Stop any still-running job for the same output folder before resuming; editing
this script does not stop an already-running IQ-TREE process. No branch-support
values are invented for the new unsupported tree.
