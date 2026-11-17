# Pipeline source notebooks

Use these source notebooks for current analyses:

| Notebook | Purpose |
| --- | --- |
| `primer_umi_seqkit_parser.ipynb` | Select best forward/reverse primer hits per raw read. |
| `sintax_output_parser.ipynb` | Parse taxonomic assignments, join primer evidence and calculate prevalence across barcode folders. |
| `evaluate_tagtracking_taxid_reduced.ipynb` | Review one barcode's read/species/primer tables and plot its species counts. |
| `sankey_target_species_detection.ipynb` | Build across-barcode summaries and plots. |

These notebooks do not load a PCR-plate mapping, create PCR-plate columns, or
group observations by PCR plate. Sample identity is `source_folder`/barcode.
PCR-primer names still identify the experimental primers and are required.
`unclassified` remains available in read tables but is excluded from barcode
prevalence. New analyses write `taxon_barcode_prevalence.tsv`.

Run `python code/read_processing.py DATASET` for the complete workflow, or
`python code/run_sintax_notebooks.py DATASET --primer-hits PATH --output-dir PATH`
to report on existing inputs. The runner executes these source notebooks and
saves report copies separately; choose a new output directory for a new analysis.

## Historical executed notebooks

`executed_notebooks/` contains local analysis artifacts and is excluded from Git.
The 46 October 2026 executed reports containing obsolete PCR-plate code were
preserved byte-for-byte in
`executed_notebooks/historical_pcr_reports_20261021.tar.gz`; the adjacent
`.manifest.json` records their original paths and SHA256 hashes. Their original
outputs remain inside that archive. They are historical records, not current
pipeline entry points. Use the source notebooks above for new work.

Postprocessing uses one shared grouping/clustering helper for local individuals
and reference taxa. Before local clustering it retains the largest extracted
taxonomic-label count per barcode/forward primer (summed across reverse primers;
all tied highest-count labels are retained). The extraction TSV records reverse primers and retention
decisions. Local reads use 97% VSEARCH identity; references use 98% and are grouped
by header taxon. See [the Postprocessing guide](../docs/post_processing.md) for
the biological assumptions and output files.

BOLD reference extraction retains named species identifications only; genus-only,
missing and provisional (`sp.`, `cf.`, `aff.`) identifications are excluded.
