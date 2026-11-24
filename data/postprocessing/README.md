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

For extracted read abundances, open `extraction_by_barcode_primer_taxon.tsv`.
For the phylogeny, open `P.figure.pdf` or `P.figure.png`. To inspect aligned bases,
open `P.aln.fasta` and `P.clipkit.aln.fasta`. For the tree model and inference
report, open `P.iqtree`. For missing sequences, inspect `translation_lengths.tsv`
and `P.kcalign_selection.tsv` when those stages were used. For execution errors,
start with `postprocessing.log` and the relevant stage's log.

Postprocessing reads the saved SINTAX classification TSVs, cleaned reads and
completed analysis; it does not need the original classifier database. Keep those
classification outputs even if the source database is not present.

## Local read extraction and clustering

| File or directory | Producer | Contents and what you can learn |
| --- | --- | --- |
| `cleaned_TAXON_DATE.fasta` | `export_sintax_taxa.py` / `export_sintax_fasta.py`, launched by Postprocessing | Actual selected cleaned query sequences before clustering. Headers retain barcode/folder, forward primer, reverse primer, lowest trusted taxon, contamination flag and cleaned read ID. These are cleaned reads, not raw reads or cluster consensuses. |
| `extraction_by_barcode_primer_taxon.tsv` | Postprocessing, from the preceding FASTA | All extracted read counts by barcode, forward primer, reverse primer and taxonomic label, plus `kept_for_module2` recording the dominant-label decision. |
| `extraction_counts.tsv` | Postprocessing | Total extracted reads for each source folder and a TOTAL row. Its total should equal the sum of `read_count` in the grouped extraction table. |
| `extraction.log` | Exporter, captured by Postprocessing | Selection/extraction messages, per-folder counts and export errors. |
| `cleaned_TAXON_DATE.dominant.fasta` | Postprocessing | Original local reads retained by dominant-taxon selection before clustering. |
| `cleaned_TAXON_DATE.vsearch97.fasta` | VSEARCH, with headers restored by Postprocessing | Local cluster consensus sequences. The header identifies the cluster seed; the sequence is a consensus and its count is not the original number of reads. |
| `vsearch_local_groups.tsv` | Postprocessing | One row per barcode/folder + forward primer + retained trusted-taxon partition; input and consensus counts and the recorded identity setting. |
| `vsearch_local_memberships.tsv` | VSEARCH UC records combined by Postprocessing | Cluster membership with local group IDs. Trace reads into the clusters used to construct consensuses. |
| `vsearch_local_groups/` | Postprocessing and VSEARCH | Per-partition input FASTAs, consensus FASTAs, UC membership files and tool logs. These explain the aggregate clustering tables. |
| `clustering_counts.tsv` | Postprocessing | Selected local and length-filtered reference input/consensus totals and the recorded identity setting; compare the two independently clustered sources. |

### Interpreting the grouped extraction table

Local extraction first keeps all requested reads in `cleaned_TAXON_DATE.fasta`.
The audit table `extraction_by_barcode_primer_taxon.tsv` has columns `barcode`,
`forward_primer`, `reverse_primer`, `taxonomic_identification`, `read_count` and
`kept_for_module2` (`true`/`false`). Counts are split by reverse primer, but selection
sums each taxonomic label across reverse primers within each barcode + forward
primer. Every label tied for the largest total survives. All reads with any of those
labels are retained regardless
of reverse primer or contamination flag; other labels are excluded from downstream
analyses. Selection uses only the requested taxon's extracted reads, not reads
outside the extraction. `none`, `unknown` and combined primer annotations remain
exact distinct forward-primer values, including in `unclassified` when requested.

The assumption is one insect individual per barcode + forward primer, whose DNA
should dominate. The most abundant label need not exceed 50% and is not proof of
identity or absence of contamination. For example, barcode15/bc1009 with Tipula
tipula (37), Tipula spectra (10) and Tipula nigris (2) retains the 37 Tipula tipula
reads. `cleaned_TAXON_DATE.dominant.fasta` contains the retained original reads;
the original extraction and every audit-table row remain available.

Counts represent cleaned queries, including separate adapter-split fragments,
not distinct raw parents or final tree tips. Reverse-primer annotations in a cluster
header describe its seed; contributing reads may have different reverse primers.
A compatible resume regenerates the table when the original extraction exists.

Retained local reads still cluster at 97% within each barcode/forward-primer
individual: a broad genus or higher-rank assignment can hide multiple species,
so all reads with that label must not automatically form one consensus. BOLD
extraction retains species-identified records only; BOLD/NCBI clustering is restricted to the same header
taxon label to preserve species positions in the tree. Local and reference
clusters are always separate.

Current local runs use `--id 0.97` and `.vsearch97.fasta`; reference runs use
`--id 0.98` and `.vsearch98.fasta`. Consult recorded commands for historical runs.

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
| `R.fasta` | Postprocessing | Reference sequences retained by the title policy. New runs append a safe taxon suffix to reference IDs so taxonomy carries into alignment/tree headers. |
| `R.length500-2000.fasta` | SeqKit | Reference sequences with inclusive lengths 500–2,000 nucleotides. |
| `R.vsearch98.fasta` | VSEARCH and Postprocessing | Reference cluster consensuses, clustered independently of local reads. Seed headers identify provenance. |
| `vsearch_reference_groups.tsv` | Postprocessing | One group per exact header taxon label; input and consensus counts at 98% identity. Different species labels cannot share clusters. |
| `vsearch_reference_memberships.tsv` | Postprocessing and VSEARCH | Group ID plus ten UC fields; trace reference contributors into their taxon-specific clusters. Cluster numbers are local to each group. |
| `vsearch_reference_groups/` | Postprocessing and VSEARCH | Per-taxon inputs, consensuses, standard UC files and logs. Older runs may instead have `R.vsearch98.uc` and a single reference log. |

BOLD and NCBI retrieval files are alternatives, so a BOLD run normally lacks the
NCBI XML files. Older runs may retain original, shorter reference identifiers;
current plotting supports both identifier formats.

## Combined sequences, translation screening and alignments

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `combined_TAXON_DATE.fasta` | Postprocessing | Local consensuses and reference consensuses combined after separate clustering. This is the input to translation screening or unfiltered selection. |
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
| `clipkit.stdout.log` | ClipKIT output captured by Postprocessing | Trimming mode, column-removal report and diagnostics. Consult the command/log for the actual mode used by that version. |

## IQ-TREE inference files

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `P.treefile` | IQ-TREE | Maximum-likelihood Newick tree with estimated branch lengths; existing runs may also contain support labels. The plotted midpoint root is written separately. |
| `P.iqtree` | IQ-TREE | Main scientific report: model selection, alignment statistics, likelihoods, inference settings and support interpretation when applicable. |
| `P.log` | IQ-TREE | Tree-search progress, warnings, timing and other inference diagnostics. |
| `iqtree.stdout.log` | IQ-TREE output captured by Postprocessing | Console output for the IQ-TREE invocation. |
| `P.bionj` | IQ-TREE | BIONJ starting/distance tree used during inference. |
| `P.mldist` | IQ-TREE | Pairwise maximum-likelihood distances. This matrix can be very large; it is not a read-abundance table. |
| `P.model.gz` | IQ-TREE | Compressed model-selection results/cache. |
| `P.ckp.gz` | IQ-TREE | Compressed checkpoint for resuming IQ-TREE itself. It is not a sequence database. |
| `P.uniqueseq.phy` | IQ-TREE | Internal PHYLIP alignment of unique sequences used by that IQ-TREE run. Inspect the IQ-TREE report for how identical input sequences were handled. |
| `P.contree`, `P.boottrees`, `P.ufboot`, `P.splits.nex` | IQ-TREE, when bootstrap/support analysis was requested | Consensus tree, bootstrap replicate trees or split/support information. Current default runs disable bootstrap, so these may exist only in older/alternative runs. |

## Rooted tree figures and labels

| File | Producer | Contents and what you can learn |
| --- | --- | --- |
| `P.figure.pdf` | Baltic + Matplotlib via `plot_phylogeny.py` | Vector phylogeny for reading, printing and sharing. |
| `P.figure.png` | Same | Raster phylogeny for quick viewing or presentations. Large trees produce tall images. |
| `P.figure.svg` | Same | Editable vector phylogeny. |
| `P.figure.midpoint.nwk` | Baltic via `plot_phylogeny.py` | Midpoint-rooted tree used for the figure, with branch lengths and preserved existing support labels. Midpoint rooting does not establish the biological root. |
| `P.figure.tips.tsv` | `plot_phylogeny.py` | Maps original tree tip IDs to source, taxon/organism, accession or BOLD process ID, local read/folder and display label. This is a tip map, not an abundance table. |
| `P.figure.plot.json` | `plot_phylogeny.py` | Plot provenance: Baltic version, input tree, rooting, number of tips, support-label handling and scale-bar units. |
| `tree_figure.stdout.log` | Plotting output captured by Postprocessing | Which tree/figure files were used or written, plus plotting diagnostics. |

## Settings, completion and resume records

| File or directory | Producer | Contents and what you can learn |
| --- | --- | --- |
| `parameters.json` | Postprocessing | Original options, data/reference choices, dependency paths and policy versions. Use this to reproduce the run. |
| `resume_parameters.json` | Postprocessing | Options from the latest resume attempt. Compare with the original settings. |
| `postprocessing.log` | Postprocessing and captured tool stderr | Executed commands, stage decisions, counts, timing/progress messages and errors. Best starting point for an execution audit. |
| `checkpoints.json` | Postprocessing | Recorded output hashes for reusable completed stages. |
| `invalidated_stages.json` | Postprocessing | Stages that must be rebuilt following a change or interrupted rebuild; an empty list means none are marked invalid. |
| `alignment.complete.json` | Postprocessing | Chosen aligner and hashes validating the alignment and retained input FASTA. |
| `clipkit.complete.json` | Postprocessing | Trimming settings/input hash used to decide whether trimming can be reused. |
| `phylogeny.complete.json` | Postprocessing | IQ-TREE model/bootstrap settings and trimmed-alignment input hash. |
| `postprocessing.complete` | Postprocessing | Success marker written after the full run finishes. |
| `postprocessing.failed` | Postprocessing, on failure | Last pipeline error. Removed after a successful retry. |
| `resume_archive/` | Postprocessing | Earlier/stale outputs moved aside before rebuilding stages. Subfolders identify when and which stage was archived. |
| `reference_header_backups/` | One-off header update, when present | Backups and mappings from a manual reference-header annotation update in an older run. Not a routine pipeline output. |
| Hidden temporary files/directories or `*.partial` | Postprocessing/tool staging | In-progress files. Their existence does not prove that a stage completed. |

Keep settings, metadata and completion records with their run. They support
resume validation and traceability; do not edit them to make an incomplete stage
appear complete.

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
[complete review instructions](../../docs/post_processing.md#review-the-trimmed-alignment-before-interpreting-the-tree).

For a field-by-field explanation of local read headers and tracing a sequence
through counts, clusters and alignments, see [Section C of the guide](../../docs/post_processing.md#common-results).
After reviewing and editing the alignment, estimate the final tree with bootstrap
support; the initial automatic pipeline run remains without bootstrap. The guide
provides a separate final IQ-TREE example and links to the
[Baltic gallery](https://phylo-baltic.github.io/baltic-gallery/) for Python figures.
