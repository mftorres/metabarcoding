"""Summary tables and figures written at the end of a Read processing analysis."""
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


def join_lists(values):
    return ';'.join(sorted({name for value in values if isinstance(value, str) for name in value.split(';') if name}))


def species_plot(species, folder, path):
    """Top 25 species with confident reads in one barcode, coloured by prevalence flag."""
    data = species[species.read_count.gt(0)].nlargest(25, 'read_count').sort_values('read_count')
    fig, ax = plt.subplots(figsize=(11, max(4, .3 * len(data) + 1)))
    ax.barh(data.species_name, data.read_count, color=data.contamination.map(FLAG_COLOURS))
    ax.set(xlabel='Cleaned reads with species confidence >=0.8', title=f'{folder}: species identification')
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_report(analysis_dir):
    analysis_dir = Path(analysis_dir)
    params = json.loads((analysis_dir/'analysis_parameters.json').read_text())
    species = pd.read_csv(analysis_dir/'all_species_summary.tsv', sep='\t',
                          dtype={'forward_primers_detected': str, 'reverse_primers_detected': str})
    qc = pd.read_csv(analysis_dir/'barcode_analysis_summary.tsv', sep='\t', keep_default_na=False)
    tags = pd.read_csv(analysis_dir/'primer_tag_summary.tsv', sep='\t')
    for folder in params['analysed_folders']:
        figures = analysis_dir/folder/'figures'
        figures.mkdir(parents=True, exist_ok=True)
        species_plot(species[species.source_folder.eq(folder)], folder, figures/f'{folder}_species_read_counts.png')
    # Confident species read counts: one row per barcode, one column per species.
    counts = species.pivot_table(index='source_folder', columns='species', values='read_count',
                                 aggfunc='sum', fill_value=0)
    counts = counts.join(qc.set_index('source_folder')[['forward_primers_detected', 'reverse_primers_detected']])
    counts.to_csv(analysis_dir/'species_by_folder.tsv', sep='\t', mode='x')
    flags = (species[species.read_count.gt(0)].groupby(['source_folder', 'contamination'])
             .agg(species_count=('species', 'nunique'), read_count=('read_count', 'sum'),
                  forward_primers_detected=('forward_primers_detected', join_lists),
                  reverse_primers_detected=('reverse_primers_detected', join_lists)).reset_index())
    flags.to_csv(analysis_dir/'contamination_summary_by_folder.tsv', sep='\t', index=False, mode='x')
    fig, ax = plt.subplots(figsize=(12, 5))
    qc.plot.bar(x='source_folder', y=['assigned_query_count', 'species_query_count'], ax=ax)
    ax.set(ylabel='Cleaned reads', xlabel='Sample folder', title='SINTAX confident assignments by folder')
    fig.tight_layout()
    fig.savefig(analysis_dir/'barcode_assignment_counts.png', dpi=160)
    plt.close(fig)
    # Tag identification outcome per folder, as a share of raw reads.
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, direction in zip(axes, ('forward', 'reverse')):
        table = (tags[tags.direction.eq(direction)].pivot_table(index='source_folder', columns='status',
                 values='raw_reads', aggfunc='sum', fill_value=0))
        table = table.reindex(columns=[s for s in STATUS_COLOURS if s in table.columns])
        (table.div(table.sum(axis=1), axis=0) * 100).plot.bar(
            stacked=True, ax=ax, color=[STATUS_COLOURS[s] for s in table.columns], legend=False)
        ax.set(title=f'{direction} tag', xlabel='Sample folder', ylabel='% of raw reads')
    axes[1].legend(handles=[Patch(color=c, label=s) for s, c in STATUS_COLOURS.items()], title='Tag status',
                   bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    fig.savefig(analysis_dir/'primer_tag_status.png', dpi=160)
    plt.close(fig)
    print(f'Report figures and summary tables written to {analysis_dir}', flush=True)
