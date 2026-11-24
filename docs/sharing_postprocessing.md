# Sharing Module 2 inputs with a collaborator

A collaborator can run Module 2 for **any taxon** in your samples without running
Module 1, if you send them the cleaned reads and the analysis tables. Raw reads,
SINTAX outputs and the MIDORI database are not needed.

## What to send

From one completed Module 1 analysis (`ANALYSIS` below):

| Files | Why |
| --- | --- |
| `ANALYSIS/analysis.complete`, `analysis_parameters.json`, `input_manifest.tsv`, `taxon_barcode_prevalence.tsv` | Mark the analysis as complete; list samples, settings and input files. |
| `ANALYSIS/barcodeNN/barcodeNN_reads.tsv`, one per sample | Taxonomy, tags and contamination flag of every read; used to select and label reads. |
| `data/DATASET/barcodeNN/filtered/asinput.fastq.gz`, one per sample | The cleaned reads themselves. |

The collaborator also needs:

- Linux or Windows with WSL2, conda, and the code (`git clone`, see [setup](setup.md));
- only the `taxonomy` conda environment (the `cutadapt` environment is for Module 1);
- the reference database they want to use (BOLD file, or internet for `--download`);
- for kc-align, a COI guide (you can add your files from `data/postprocessing/coi_refs/`);
- optionally UGENE and FigTree to review the alignment and tree.

## Make the archive

Run from the project folder. The first command writes the list of files; the second
packs them. Nothing is analysed again.

```bash
ANALYSIS=data/read_processing/DATASET_run01/analysis
python - "$ANALYSIS" > share.list <<'PY'
import json, sys
from pathlib import Path
analysis = Path(sys.argv[1])
params = json.loads((analysis/'analysis_parameters.json').read_text())
for name in ('analysis.complete', 'analysis_parameters.json', 'input_manifest.tsv', 'taxon_barcode_prevalence.tsv'):
    print(analysis/name)
for folder, files in params['input_files'].items():
    print(analysis/folder/f'{folder}_reads.tsv')
    print(files['clean_reads']['path'])
PY
tar -czf postprocessing-inputs.tar.gz --verbatim-files-from --files-from=share.list
```

Every path in `share.list` must start with `data/`; if one starts with `/`, the
analysis was made outside the project folder and the collaborator would need the
same absolute path. Add `data/postprocessing/coi_refs/*.fasta` to the `tar` command
to include your COI guides. The reads are already compressed, so expect an archive
of roughly the same size as the files.

If a file-size limit applies, split the archive and send all parts:

```bash
split -b 900M -d -a 3 postprocessing-inputs.tar.gz postprocessing-inputs.tar.gz.part-
```

## What to tell the collaborator

1. Follow the [setup guide](setup.md) (software and project).
2. Extract the archive in the project folder:
   [Module 2, A.1](post_processing.md#a1-a-completed-module-1-analysis).
3. Prepare the references ([BOLD](../data/bold/README.md)) and, for kc-align, a
   [COI guide](coi_reference.md).
4. Run Module 2 with `--analysis-dir` set to the `ANALYSIS` path you used above:
   [Module 2, B](post_processing.md#b-run-module-2).
