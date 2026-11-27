"""Species-by-folder table and figures written at the end of a Module 1 analysis."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

FLAG_COLOURS = {'contamination': '#c55a45', 'nocontamination': '#3486a0'}
STATUS_COLOURS = {'assigned': '#3486a0', 'ambiguous': '#e0a43a', 'too_many_errors': '#c55a45',
                  'truncated': '#9b7fb8', 'primer_not_found': '#b7c1c5'}


def species_plot(species, folder, path):
    """The 25 species with most reads in one folder, coloured by the contamination flag."""
    data = species.nlargest(25, 'read_count').sort_values('read_count')
    fig, ax = plt.subplots(figsize=(11, max(4, 0.3 * len(data) + 1)))
    ax.barh(data.species_name, data.read_count, color=data.contamination.map(FLAG_COLOURS))
    ax.set(xlabel='Cleaned reads in clusters named to species', title=f'{folder}: species identification')
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_report(analysis_dir):
    analysis_dir = Path(analysis_dir)
    params = json.loads((analysis_dir/'analysis_parameters.json').read_text())
    species = pd.read_csv(analysis_dir/'all_species_summary.tsv', sep='\t', dtype={'forward_tags': str})
    stats = pd.read_csv(analysis_dir/'barcode_analysis_summary.tsv', sep='\t')
    tags = pd.read_csv(analysis_dir/'primer_tag_summary.tsv', sep='\t')

    # Species read counts: one row per folder, one column per species.
    counts = species.pivot_table(index='source_folder', columns='species', values='read_count',
                                 aggfunc='sum', fill_value=0)
    counts.to_csv(analysis_dir/'species_by_folder.tsv', sep='\t', mode='x')

    for folder in params['analysed_folders']:
        species_plot(species[species.source_folder.eq(folder)], folder,
                     analysis_dir/folder/f'{folder}_species_read_counts.png')

    # Reads identified per folder.
    fig, ax = plt.subplots(figsize=(12, 5))
    stats.plot.bar(x='source_folder', y=['cleaned_reads', 'reads_assigned', 'reads_to_species'], ax=ax)
    ax.set(ylabel='Cleaned reads', xlabel='Sample folder',
           title=f'Reads identified with confidence >= {params["confidence_cutoff"]}')
    fig.tight_layout()
    fig.savefig(analysis_dir/'barcode_assignment_counts.png', dpi=160)
    plt.close(fig)

    # Tag identification outcome per folder, as a share of raw reads.
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, direction in zip(axes, ('forward', 'reverse')):
        table = tags[tags.direction.eq(direction)].pivot_table(index='source_folder', columns='status',
                                                                 values='raw_reads', aggfunc='sum', fill_value=0)
        table = table.reindex(columns=[status for status in STATUS_COLOURS if status in table.columns])
        (table.div(table.sum(axis=1), axis=0) * 100).plot.bar(
            stacked=True, ax=ax, color=[STATUS_COLOURS[status] for status in table.columns], legend=False)
        ax.set(title=f'{direction} tag', xlabel='Sample folder', ylabel='% of raw reads')
    axes[1].legend(handles=[Patch(color=colour, label=status) for status, colour in STATUS_COLOURS.items()],
                   title='Tag status', bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    fig.savefig(analysis_dir/'primer_tag_status.png', dpi=160)
    plt.close(fig)
    print(f'Report figures and species_by_folder.tsv written to {analysis_dir}', flush=True)
