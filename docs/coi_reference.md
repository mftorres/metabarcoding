# Choosing a COI reference for kc-align

**For COI, kc-align is the preferred alignment choice in this project when a
trustworthy coding reference is available.** MAFFT remains the software default:
leaving off `--kcalign` still runs MAFFT. For noncoding DNA, use MAFFT rather than
kc-align. The coding guide is a separate input from the BOLD/NCBI comparison
sequences downloaded for the tree; choosing `--download` does not prepare the
coding guide automatically.

## 1. Choose a trustworthy source

Use an identified specimen from the same species, if available, or a closely
related species in the same genus or family that covers the same COI region.
Prefer a curated or published record with specimen/voucher information and
supporting evidence that it is mitochondrial COI. Do not automatically use the
first search result, an unverified sequence, or a local consensus you are trying
to identify. A file already in `coi_refs/` is not proof of suitability for your taxon.

**From NCBI:** search [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) for
your group and COI, COX1 or cytochrome c oxidase subunit I. Inspect the GenBank
record, its mitochondrial source, publication, specimen and **CDS** annotation
(the protein-coding DNA feature). Check its strand, exact coordinates,
`codon_start`, translation and genetic code. Download the nucleotide sequence of
that COI feature or a verified target-region extract, not the protein FASTA or an
entire mitochondrial genome containing other genes. An annotated mitochondrial
genome can be a source for extracting COI even though whole genomes are normally
excluded from Module 2's comparison-reference download. See
[NCBI's explanation of GenBank features](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/).

**From BOLD:** select a named-species COI record with credible specimen
identification, voucher/institution information and, where available, trace data
and a supporting publication or linked GenBank record. Check the marker and
sequence quality, and download the nucleotide FASTA. The bulk database header
alone does not establish the exact coding boundaries or rule out a NUMT. If the
BOLD barcode lacks coding annotation, align it to a trustworthy annotated
mitochondrial COI sequence to establish its frame and precise target interval,
or ask the maintainer for a prepared guide.

## 2. Exclude suspected NUMTs

A **NUMT** is a copy of mitochondrial DNA present in the nuclear genome. The guide
should be genuine mitochondrial COI, not a NUMT, pseudogene, contaminant or chimera.
Check the source annotation and evidence of mitochondrial origin. Compare its
protein translation and nucleotide alignment with independently trusted COI
records. Reject suspected NUMTs, unexpected frameshifts, internal stops or
unexplained insertions/deletions.

**No internal stops is necessary but not sufficient.** Some NUMTs retain an
apparently intact reading frame. A database name or a successful pipeline check
cannot certify that the guide is NUMT-free. If the evidence is uncertain, choose a
better documented record rather than calling it verified. This limitation is
supported by [research on insect COI NUMTs](https://pmc.ncbi.nlm.nih.gov/articles/PMC10249859/).

## 3. Prepare the exact COI coding region, without primers

**Use the exact expected COI open reading frame (ORF) or the precisely defined
coding interval being analysed.** Establish its boundaries from annotation and
known codon positions, not by guessing a length. There is no single universal COI
length for all taxa and all assays.

- If your analysis uses the **complete COI ORF**, use the complete annotated COI
  coding sequence of the expected length for that record/taxon, in the correct
  orientation. A short COI-5P barcode is not a complete COI ORF.
- If your reads cover only a **COI barcode fragment**, use a verified, in-frame
  extract corresponding to that same region. Module 1 cuts the primers off, so for
  this project's primers the reads span the ~658 nt **between** the forward primer
  (`…TCAACMAATCATAAAGATATTGG`) and the reverse primer (`TAAACTTCWGGRTGWCCAAAAAATCA`).
  Use a guide covering that region. It must start and finish on complete
  codon boundaries, with its exact expected coding length documented. Do not use
  an unrelated fragment or assume a fixed barcode length establishes the frame.
- Remove laboratory primers, primer tails, tags, adapters and flanking non-COI
  sequence. Use the known primer positions and coding annotation to define the
  interval; do not remove genuine COI bases solely because a short motif resembles
  a primer. The guide must contain only the intended biological coding region.
- Save **one nucleotide sequence** in a plain FASTA file: one `>` header followed
  by ungapped A/C/G/T DNA, without `N` or other ambiguity codes. Keep it in the
  coding orientation. Its length must be divisible by three, and translation with
  the correct genetic code must contain no internal stop codons.
- Do not add or delete arbitrary bases to make the length divisible by three.
  Resolve partial codons or incomplete mitochondrial stop annotations using the
  source annotation and a validated preparation, not automatic padding. A genuine
  terminal stop must agree with the annotated complete CDS; a barcode fragment
  does not need a start or terminal stop codon.

For the usual insect mitochondrial COI analysis, use genetic code **5**. Confirm
the correct code for other organisms. If you cannot verify the interval or its
reading frame, obtain a prepared reference from the maintainer before proceeding.

## 4. Save and check the guide

From the project root, create its folder:

```bash
mkdir -p data/postprocessing/coi_refs
```

Save the prepared FASTA as `data/postprocessing/coi_refs/COI_reference.fasta`, or
choose a descriptive name and use that path in every command. Renaming alone does
not remove primers, select the right region, or fix a reading frame. Keep the
original downloaded record separately. Beside the guide, record the accession
with version or BOLD process ID, download date, source link, species, coordinates,
strand, expected coding length, genetic code, primer removal and evidence used to
exclude NUMTs. These notes are important when someone else repeats the analysis.

You can run the pipeline's basic checks before a long analysis:

```bash
python -c "import sys; sys.path.insert(0, 'code'); from kcalign_alignment import validate_reference; name, sequence = validate_reference('data/postprocessing/coi_refs/COI_reference.fasta', table=5); print(name); print('Length:', len(sequence), 'nt'); print('Basic checks passed; verify source, exact COI boundaries and absence of NUMTs separately.')"
```

This checks the number of sequences, letters, length divisibility and internal
stops. It does **not** verify species identity, exact ORF boundaries, primer
absence or mitochondrial origin. Compare the reported length with the expected
length established from the source record, and check the alignment after the run.

## 5. Run Module 2 with kc-align

After activating the environment and preparing the guide, run this example from
the project root. Replace the taxon, barcode and analysis path as needed:

```bash
python code/postprocess_sintax.py --bold \
  --bold-db data/bold/BOLD_filtered.fasta.gz \
  --taxon Tipulidae --folders barcode15 \
  --analysis-dir data/read_processing/DATASET_run01/analysis \
  --output-dir data/postprocessing/Tipulidae_barcode15_kcalign_run01 \
  --kcalign --kcalign-reference data/postprocessing/coi_refs/COI_reference.fasta \
  --translation-table 5 --threads 6
```

For NCBI comparison sequences, replace `--bold` with `--download` and remove
`--bold-db ...`; keep the kc-align options and the independently prepared guide.
Use **1 thread or a multiple of three**, such as 6. Leave the protein screen on
for the usual COI run. kc-align may reject problematic coding sequences; check
`*.kcalign_selection.tsv` and the log. At least four aligned consensus sequences
must remain for a tree. The coding guide is not automatically added as a tree tip.

If you received a shared archive instead of running Module 1, use the analysis path
the sender gave you as `--analysis-dir` (see
[Module 2, A.1](post_processing.md#a1-a-completed-module-1-analysis)).

## Review and final record selection

For COI comparison references, prefer BOLD as the starting source; NCBI remains
a useful alternative and a source of well-annotated coding guides. Neither source
guarantees record quality. You must inspect the ClipKIT-trimmed alignment in
[UGENE](https://ugene.net/download-all.html) before interpreting the tree. Use
[Peartree](https://peartree.live/) to quickly view `*.figure.midpoint.nwk`.
Final decisions about local reads/consensuses, BOLD species and individual
reference records are at the user’s scientific discretion. Document exclusions
and re-estimate the tree after changing its alignment or sequence set. See the
[full review instructions](post_processing.md#c-review-the-results).
