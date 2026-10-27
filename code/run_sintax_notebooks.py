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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--primer-hits',type=Path,default=Path('data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_20260731.txt'))
    parser.add_argument('--folders',nargs='+',default=[f'barcode{i:02d}' for i in range(1,25)])
    parser.add_argument('--resume',action='store_true',help='Reuse a completed parsing stage and execute remaining reports')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    output=args.output_dir.resolve()
    if output.exists() and not args.resume: raise FileExistsError(output)
    if args.resume and not (output/'analysis.complete').exists(): raise ValueError('Cannot resume incomplete parsing')
    executed=root/'code/executed_notebooks'/('sintax_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    executed.mkdir(parents=True)
    os.environ['SINTAX_ANALYSIS_ROOT']=str(output)
    os.environ['SINTAX_FOLDERS']=','.join(args.folders)
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
