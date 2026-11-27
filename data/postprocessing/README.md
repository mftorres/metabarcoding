# Postprocessing run folders and output files

Each run folder, such as `Ceratopogonidae_bold_kcalign_run01/`, contains one
postprocessing analysis started by `code/postprocess_sintax.py`. This guide
explains the files in that example and their equivalents in other runs.
Run settings determine which files appear: BOLD versus NCBI, MAFFT versus
kc-align, translation screening versus `--skip-translation`, and IQ-TREE options.
Database downloads are described in [the BOLD README](../bold/README.md) and the
[main Postprocessing guide](../../docs/post_processing.md).

In the tables below:

- `TAXON_DATE` stands for e.g. `Ceratopogonidae_20261030`.
- `R` stands for `BOLD_TAXON_DATE` or `COI_TAXON_DATE`.
- `P` stands for `combined_TAXON_nonumt_DATE` (translation-screened) or
  `combined_TAXON_selected_DATE` (translation skipped).
- `.tsv` files are tab-delimited tables with column headers. FASTA files contain
  `>sequence_identifier` lines followed by sequences; `.aln.fasta` contains gaps.

## Which files to open first

For the reads and clusters of each specimen, open `extraction_by_specimen_taxon.tsv`.
For the phylogeny, open `P.figure.pdf` or `P.figure.png`. To inspect aligned bases,
open `P.aln.fasta` and `P.clipkit.aln.fasta`. For the tree model and inference
report, open `P.iqtree`. For missing sequences, inspect `translation_lengths.tsv`
and `P.kcalign_selection.tsv` when those stages were used. For execution errors,
start with `postprocessing.log` and the relevant stage's log.

Postprocessing reads the cluster tables and cluster sequences of the completed
Module 1 analysis; it needs neither the SINTAX outputs nor the MIDORI database.

## Local clusters

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `clusters_TAXON_DATE.fasta` | `export_clusters.py`, launched by Postprocessing | Every Module 1 cluster (consensus sequence or single read) accepted as the taxon in the chosen folders. The names give folder, forward tag, identification, contamination flag, read count and cluster ID ([details](../../docs/post_processing.md#reading-the-local-sequence-names)). |
| `extraction.log` | Exporter, captured by Postprocessing | Clusters exported per folder, and export errors. |
| `extraction_by_specimen_taxon.tsv` | Postprocessing | One row per specimen (`source_folder` + `forward_tag`) and identification (`taxon`): `cluster_count`, `read_count`, and `kept_for_tree` (`true` for the identification with most reads in the specimen). |
| `clusters_TAXON_DATE.main_taxon.fasta` | Postprocessing | The clusters kept: those of the main identification of each specimen. These go on to the protein check, alignment and tree. |

### How the main identification is chosen

Each PCR well (barcode folder + forward tag) holds one arthropod, whose DNA should
give most reads. Within each specimen, the identifications of the exported clusters
are compared by their total reads and the largest is kept, with all its clusters;
tied identifications are all kept. For example, barcode15/bc1009 with *Tipula
tipula* (37 reads), *T. spectra* (10) and *T. nigris* (2) keeps the clusters of
*T. tipula*. The largest identification need not exceed 50% and is not proof of
identity or of the absence of contamination. Only clusters of the requested taxon
are compared, not the specimen's other reads (see `specimen_taxa.tsv` in the Module 1
analysis for those). `none` means no forward tag was assigned; such a group can mix
specimens.

A specimen can keep several clusters of its main identification: reads less than 97%
identical to each other, often single reads with many sequencing errors, form
separate clusters. Each becomes a separate tip.

## External references: BOLD or NCBI

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `reference_query.txt` | Postprocessing | BOLD header-match text or the NCBI Entrez query defining reference retrieval. |
| `references_bold.fasta` | `bold_references.py` | Taxon-matching records selected from the local BOLD database, before reference filtering and clustering. |
| `bold_reference_metadata.tsv` | `bold_references.py` | Maps generated BOLD sequence IDs to process IDs, marker, available taxon, original header and source database. Use it to recover provenance or look up a reference in BOLD. |
| `taxonomy_search.xml` | NCBI Entrez Direct `esearch` | Taxonomy search response when an NCBI taxon name was supplied. |
| `taxonomy.xml` | NCBI Entrez Direct `efetch` | Resolved taxon name and taxid; confirms what the NCBI query means. |
| `coi_search.xml` | NCBI Entrez Direct `esearch` | NCBI nucleotide search response used to fetch matching COI references. |
| `references_downloaded.fasta` | NCBI Entrez Direct `efetch` | Downloaded nucleotide reference records before filtering. |
| `reference_records.gb.xml` | NCBI Entrez Direct `efetch` | GenBank organism/accession metadata for reference labels and plotting. Entrez may concatenate multiple XML batches. |
| `reference_accessions.txt` | BOLD extractor or NCBI Entrez Direct | IDs of references initially selected/downloaded; BOLD uses its generated sequence IDs. |
| `reference_filter.tsv` | Postprocessing | Reference descriptions, lengths, retain/exclude decisions and reasons, including the mitogenome-title policy. |
| `reference_accessions_retained.txt` | Postprocessing | IDs surviving that reference policy, before the later length filter. |
| `R.fasta` | Postprocessing | Reference sequences retained by the title policy, with a `__taxon__Species_name` suffix on each ID so the species carries into alignment and tree. |
| `R.length500-2000.fasta` | SeqKit | Reference sequences with inclusive lengths 500–2,000 nucleotides. |
| `R.vsearch98.fasta` | VSEARCH and Postprocessing | Reference cluster consensuses, clustered independently of local reads. Seed headers identify provenance. |
| `vsearch_reference_groups.tsv` | Postprocessing | One group per exact header taxon label; input and consensus counts at 98% identity. Different species labels cannot share clusters. |
| `vsearch_reference_memberships.tsv` | Postprocessing and VSEARCH | Group ID plus ten UC fields; trace reference contributors into their taxon-specific clusters. Cluster numbers are local to each group. |
| `vsearch_reference_groups/` | Postprocessing and VSEARCH | Per-taxon inputs, consensuses, standard UC files and logs. |

BOLD and NCBI retrieval files are alternatives, so a BOLD run lacks the NCBI XML files.
`bold_reference_metadata.tsv` has a `marker` column: the BOLD file also contains
markers other than COI-5P (see [the BOLD README](../bold/README.md#other-markers-in-the-bold-file)).

## Combined sequences, translation screening and alignments

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `combined_TAXON_DATE.fasta` | Postprocessing | Local clusters and reference consensus sequences together. This is the input to translation screening or unfiltered selection. |
| `combined_TAXON_DATE.six_frames.faa` | SeqKit `translate` | Six translated reading frames per input representative under the chosen genetic code. Investigate coding potential and rejected sequences here. |
| `translation_lengths.tsv` | Postprocessing | Translation-screening evidence and keep/reject decisions for each representative; use the headers to inspect the recorded frame/length fields. |
| `P.fasta` | Postprocessing | Representatives retained for alignment: `nonumt` when screened, or `selected` when translation screening was skipped. The noNUMT name describes a screening heuristic, not proof that every NUMT was excluded. |
| `P.aln.fasta` | MAFFT or kc-align, restored/validated by Postprocessing | Untrimmed nucleotide multiple-sequence alignment, with original identifiers restored. This is the base-level view before column trimming. kc-align may exclude some inputs. |
| `alignment_orientation.tsv` | Postprocessing after MAFFT | Original sequence IDs and alignment strand, including sequences reverse-complemented by MAFFT. Usually absent from kc-align runs. |
| `P.kcalign_ids.tsv` | Postprocessing kc-align adapter | Maps temporary tool IDs to original headers. |
| `P.kcalign_selection.tsv` | Postprocessing kc-align adapter | Records each representative as aligned or rejected by kc-align. Explains differences between retained FASTA and alignment membership. |
| `P.pre_align.aa.fasta` | kc-align via the adapter | Amino-acid pre-alignment retained for inspection. It is not the nucleotide alignment passed to IQ-TREE. |
| `P.kcalign.stdout.log` | kc-align output captured by Postprocessing | Coding alignment messages, screening/rejection information and tool diagnostics. |
| `P.clipkit.aln.fasta` | ClipKIT | Column-trimmed nucleotide alignment supplied to IQ-TREE. Compare with `P.aln.fasta` to inspect removed columns. |
| `clipkit.stdout.log` | ClipKIT output captured by Postprocessing | Trimming report: ClipKIT `gappy` mode removes columns that are ≥90% gaps. |

## IQ-TREE inference files

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `P.treefile` | IQ-TREE | Maximum-likelihood Newick tree with estimated branch lengths, without support values. The plotted midpoint root is written separately. |
| `P.iqtree` | IQ-TREE | Main scientific report: model selection, alignment statistics, likelihoods, inference settings and support interpretation when applicable. |
| `P.log` | IQ-TREE | Tree-search progress, warnings, timing and other inference diagnostics. |
| `iqtree.stdout.log` | IQ-TREE output captured by Postprocessing | Console output for the IQ-TREE invocation. |
| `P.bionj` | IQ-TREE | BIONJ starting/distance tree used during inference. |
| `P.mldist` | IQ-TREE | Pairwise maximum-likelihood distances. This matrix can be very large; it is not a read-abundance table. |
| `P.model.gz` | IQ-TREE | Compressed model-selection results/cache. |
| `P.ckp.gz` | IQ-TREE | Compressed checkpoint for resuming IQ-TREE itself. It is not a sequence database. |
| `P.uniqueseq.phy` | IQ-TREE | Internal PHYLIP alignment of unique sequences used by that IQ-TREE run. Inspect the IQ-TREE report for how identical input sequences were handled. |

## Rooted tree figures and labels

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `P.figure.pdf` | Baltic + Matplotlib via `plot_phylogeny.py` | Vector phylogeny for reading, printing and sharing. |
| `P.figure.png` | Same | Raster phylogeny for quick viewing or presentations. Large trees produce tall images. |
| `P.figure.svg` | Same | Editable vector phylogeny. |
| `P.figure.midpoint.nwk` | Baltic via `plot_phylogeny.py` | Midpoint-rooted tree used for the figure, with branch lengths. Midpoint rooting does not establish the biological root. |
| `P.figure.tips.tsv` | `plot_phylogeny.py` | Maps tree tip IDs to source, taxon/organism, accession or BOLD process ID, and for local clusters the cluster ID, read count, folder and forward tag. |
| `P.figure.plot.json` | `plot_phylogeny.py` | Plot provenance: Baltic version, input tree, rooting, number of tips, support-label handling and scale-bar units. |
| `tree_figure.stdout.log` | Plotting output captured by Postprocessing | Which tree/figure files were used or written, plus plotting diagnostics. |

## Settings, completion and resume records

| File or directory | Producer | Contents and what you can learn |
| --- | --- | --- |
| `parameters.json` | Postprocessing | Original options, data/reference choices and program paths. Use this to reproduce the run. |
| `resume_parameters.json` | Postprocessing | Options from the latest resume attempt. Compare with the original settings. |
| `postprocessing.log` | Postprocessing and captured tool stderr | Executed commands, stage decisions, counts, timing/progress messages and errors. Best starting point for an execution audit. |
| `checkpoints.json` | Postprocessing | For every completed step: output checksums, input sizes/times and settings. `--resume` reuses a step only when all three are unchanged. |
| `postprocessing.complete` | Postprocessing | Success marker written after the full run finishes. |
| `postprocessing.failed` | Postprocessing, on failure | Last pipeline error. Removed after a successful retry. |
| `resume_archive/` | Postprocessing | Earlier/stale outputs moved aside before rebuilding stages. Subfolders identify when and which stage was archived. |
| Hidden temporary files/directories or `*.partial` | Postprocessing/tool staging | In-progress files. Their existence does not prove that a stage completed. |

Keep settings, metadata and completion records with their run. They support
resume and traceability; do not edit them to make an incomplete step appear complete.
Runs made before 24 November 2026 clustered local reads in Module 2 (files
`cleaned_*.fasta`, `*.vsearch97.fasta`, `vsearch_local_*`) and cannot be resumed with
the current code.

## The shared `coi_refs/` folder

`coi_refs/` contains externally supplied trusted coding-reference FASTAs for
`--kcalign-reference`, such as `MT782148_cullicoides.fasta`,
`MZ608358_Sciaridae_sp.fasta` and `SATIP962-10_Nigrotipula_nigra.fasta`.
These are alignment guides, not extracted local reads or automatic outputs of
the run folder. The guide needs one ungapped, in-frame A/C/G/T coding sequence
under the chosen genetic code. The adapter removes its synthetic guide entry
from the final alignment.

BOLD reference extraction retains named species identifications only; genus-only,
missing and provisional (`sp.`, `cf.`, `aff.`) identifications are excluded.

## Review before interpreting the tree

The user must check `*.clipkit.aln.fasta` after ClipKIT, for example in
[UGENE](https://ugene.net/download-all.html). The pipeline does not pause before
IQ-TREE; its first tree is a draft pending this review. Quickly explore
`*.figure.midpoint.nwk` in [Peartree](https://peartree.live/).
The user decides which local reads/consensuses and reference species/records are
appropriate. Keep originals, document changes and re-estimate the tree after
changing the sequence set or alignment. See the
[complete review instructions](../../docs/post_processing.md#c-review-the-results).

For a field-by-field explanation of local sequence names, see
[Reading the local sequence names](../../docs/post_processing.md#reading-the-local-sequence-names).
After reviewing and editing the alignment, estimate the final tree with bootstrap
support; the initial automatic pipeline run remains without bootstrap. The guide
provides a separate final IQ-TREE example and links to the
[Baltic gallery](https://phylo-baltic.github.io/baltic-gallery/) for Python figures.
