#!/usr/bin/env python3
"""Execute SINTAX post-classification notebooks in the current Python environment."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import sys
import tempfile
import nbformat
from nbclient import NotebookClient
from describe_table_headers import describe_headers
from sample_folders import dataset_path, discover_samples


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--executed-dir',type=Path,help='Destination for executed report notebooks')
    parser.add_argument('dataset', nargs='?', type=Path, help='Dataset parent, e.g. sequencing2')
    parser.add_argument('--merged-dir', type=Path)
    parser.add_argument('--sintax-root', type=Path, help='Default: data/sintax/DATASET; pass data/sintax explicitly for legacy runs')
    parser.add_argument('--plate-map', type=Path, default=Path('saproxDNA/species_info_seq5/seq5_pcr_to_pools.tsv'))
    parser.add_argument('--primer-hits',type=Path,help='Prepared best-hit table; default: latest dated table for this dataset')
    parser.add_argument('--folders',nargs='+',help='Optional sample-name subset; default: discover barcodeNN and unclassified')
    parser.add_argument('--resume',action='store_true',help='Reuse a completed parsing stage and execute remaining reports')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    if args.dataset and args.merged_dir:
        parser.error('Use a dataset argument or --merged-dir, not both')
    if args.resume:
        saved = json.loads((args.output_dir/'analysis_parameters.json').read_text())
        if (args.dataset or args.merged_dir) and Path(saved['merged_dir']).resolve() != dataset_path(args.dataset or args.merged_dir):
            parser.error('Resume dataset differs from saved analysis')
        merged = Path(saved['merged_dir'])
        folders = args.folders or saved['requested_folders']
    else:
        if not (args.dataset or args.merged_dir):
            parser.error('Provide a dataset parent or --merged-dir for a new analysis')
        merged = dataset_path(args.dataset or args.merged_dir)
        folders = args.folders or [p.name for p in discover_samples([merged])]
    if args.primer_hits is None and not args.resume:
        candidates = sorted((root/'data/primertrack').glob(f'{merged.name}_seqkit_primer_umi_best_hits_*.txt'))
        if not candidates:
            parser.error('No prepared primer best-hit table found; run read_processing.py or supply --primer-hits')
        args.primer_hits = candidates[-1]
    args.sintax_root = args.sintax_root or (Path(saved['sintax_root']) if args.resume else root/'data/sintax'/merged.name)
    output=args.output_dir.resolve()
    if output.exists() and not args.resume: raise FileExistsError(output)
    if args.resume and not (output/'analysis.complete').exists(): raise ValueError('Cannot resume incomplete parsing')
    executed=args.executed_dir or root/'code/executed_notebooks'/('sintax_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    executed.mkdir(parents=True)
    os.environ['SINTAX_ANALYSIS_ROOT']=str(output)
    os.environ['SINTAX_FOLDERS']=','.join(folders)
    os.environ['SINTAX_MERGED_DIR']=str(merged.resolve())
    os.environ['SINTAX_ROOT']=str(args.sintax_root.resolve())
    os.environ['SINTAX_PLATE_MAP']=str(args.plate_map.resolve())
    if args.primer_hits is not None:
        os.environ['PRIMER_BEST_HITS']=str(args.primer_hits.resolve())
    os.environ['MPLBACKEND']='Agg'
    with tempfile.TemporaryDirectory(prefix='sintax_jupyter_') as temp:
        temp=Path(temp)
        os.environ['MPLCONFIGDIR']=str(temp/'matplotlib')
        spec=temp/'kernels/taxonomy-analysis';spec.mkdir(parents=True)
        (spec/'kernel.json').write_text(json.dumps({'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],
            'display_name':'taxonomy analysis','language':'python'}))
        os.environ['JUPYTER_PATH']=str(temp)+os.pathsep+os.environ.get('JUPYTER_PATH','')
        def execute(source,destination):
            print(f'Executing {source} -> {destination}',flush=True)
            notebook=nbformat.read(root/'code'/source,as_version=4)
            client=NotebookClient(notebook,timeout=-1,kernel_name='taxonomy-analysis',resources={'metadata':{'path':str(root)}})
            try: client.execute()
            finally: nbformat.write(notebook,executed/destination)
            # Echo concise stream logs for the batch parser; tables remain in notebooks.
            for cell in notebook.cells:
                for out in cell.get('outputs',[]):
                    if out.get('output_type')=='stream': print(out.get('text',''),end='',flush=True)
        if not args.resume: execute('sintax_output_parser.ipynb','sintax_output_parser_executed.ipynb')
        metadata=json.loads((output/'analysis_parameters.json').read_text())
        for folder in metadata['analysed_folders']:
            os.environ['SELECTED_BARCODE']=folder
            execute('evaluate_tagtracking_taxid_reduced.ipynb',f'{folder}_read_analysis.ipynb')
        if not (output/'species_by_folder.tsv').exists():
            execute('sankey_target_species_detection.ipynb','across_barcode_analysis.ipynb')
    describe_headers(output)
    (output/'notebooks.complete').write_text(str(executed)+'\n')
    print(f'Executed notebooks: {executed}',flush=True)
    print(f'Analysis outputs: {output}',flush=True)

if __name__=='__main__': main()
