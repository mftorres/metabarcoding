# Deadwood insect COI metabarcoding

This project turns Oxford Nanopore (ONT) sequencing reads of the insect COI barcode
gene into a **draft alignment and phylogenetic tree that a person then reviews**.
It is written for biologists: every step can be run by copying commands, and the
guides explain what each step does and what to check.

Open the [offline browser guide](docs/metabarcoding-guide.html) (double-click the
file; no internet or server needed) for a menu-driven version of these instructions
with copy buttons and a command builder.

## What the pipeline does

Each sequencing **barcode folder** (`barcode01`, `barcode02`, … as written by MinKNOW)
holds the reads of one pool of PCR products. Inside a pool, each insect specimen was
amplified with a different **tagged forward primer** (`bc1001`–`bc1096`) and one of three
**tagged reverse primers** (`bc1100rc`, `bc1107rc`, `bc1122rc`). The tag tells you which
specimen a read came from.

![Pipeline map: Module 1 read processing and Module 2 postprocessing, with the path each kind of dataset takes, its inputs, key outputs and main tools](docs/pipeline_diagram.svg)

Read the map left to right. Each coloured line is one stream of data: your reads (blue), the tags that identify specimens (yellow), the MIDORI classification database (purple), the BOLD or NCBI reference sequences (orange), the optional COI guide for kc-align (pink) and a shared archive for collaborators who only run Module 2 (grey). Where lines run side by side through one step, that step uses both; the local reads become one consensus per specimen, the references are reduced to one centroid per ≥98% species cluster, and the two join at "Combine". Squares are files you provide; grey pills are key outputs; the shaded area is Module 2.

The pipeline never makes the final scientific decisions. You decide which reads,
consensus sequences and reference species belong in the final tree.

## Quick start

1. Set up Ubuntu/WSL, Miniconda and the two software environments:
   [setup guide](docs/setup.md).
2. Download the classification database ([MIDORI](data/sintax/README.md)) and, for
   Module 2, the reference database ([BOLD](data/bold/README.md)).
3. Put your raw reads in `data/DATASET/barcodeNN/` (one folder per barcode; any
   dataset name works) and run, from the project folder:

```bash
conda activate taxonomy
# Only if MinKNOW wrote pass/fail subfolders: combine them first.
python code/merge_pass_fail.py DATASET
# Module 1: preview, then run.
python code/read_processing.py DATASET --dry-run
python code/read_processing.py DATASET --output-dir data/read_processing/DATASET_run01
# Module 2: one taxon.
python code/postprocess_sintax.py DATASET --bold --taxon Tipulidae \
  --analysis-dir data/read_processing/DATASET_run01/analysis \
  --output-dir data/postprocessing/Tipulidae_run01 \
  --kcalign --kcalign-reference data/postprocessing/coi_refs/COI_reference.fasta --threads 6
```

4. Review the trimmed alignment and the tree:
   [what to check](docs/post_processing.md#c-review-the-results).

## Documentation

| Guide | Read it when you want to… |
| --- | --- |
| [Setup](docs/setup.md) | install the software and get the project |
| [Module 1 · Read processing](docs/read_processing.md) | clean, tag and classify reads; understand the tables; export reads |
| [Module 2 · Postprocessing](docs/post_processing.md) | build consensus sequences, an alignment and a draft tree for one taxon, then review it |
| [Choosing a COI guide for kc-align](docs/coi_reference.md) | prepare the trusted reference sequence kc-align needs |
| [Sharing Module 2 inputs](docs/sharing_postprocessing.md) | send a collaborator what they need to run Module 2 without Module 1 |
| [Glossary](docs/glossary.md) | look up a term |
| [Module 2 output files](data/postprocessing/README.md) | find out what a file in a run folder contains |
| [Code overview](code/README.md) | see which script does what |

## Things to know before you start

- **Classification is the slow step.** SINTAX can take days for a large barcode
  folder. Everything before and after it takes minutes to hours.
- **Tags must be complete to be trusted.** A tag is assigned only when the read
  contains the whole tag with at most 4 errors and it clearly beats every other tag;
  otherwise the read is marked *ambiguous*, *too_many_errors* or *truncated*
  ([details](docs/read_processing.md#how-tags-are-identified)).
- **The BOLD database contains other markers** (16S, 28S, ITS, COI-3P, CYTB…). If you do
  not use kc-align with a COI guide, references of other markers can end up in the
  alignment ([details](data/bold/README.md#other-markers-in-the-bold-file)).
- **The first tree is a draft.** It has no bootstrap support and has not been checked
  by a person.

## Licence

See [LICENSE](LICENSE).
