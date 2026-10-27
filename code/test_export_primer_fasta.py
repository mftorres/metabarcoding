"""Primer exports tested with real SeqKit and isolated raw/clean/SINTAX fixtures."""
import csv
import gzip
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from export_primer_fasta import export_primer
from export_sintax_fasta import export_reads


def fastq(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt') as handle:
        for header, sequence in records:
            handle.write(f'@{header}\n{sequence}\n+\n'+ 'I'*len(sequence)+'\n')


@unittest.skipUnless(shutil.which('seqkit'), 'SeqKit unavailable')
class PrimerExportTests(unittest.TestCase):
    def fixture(self, root):
        analysis=root/'analysis';analysis.mkdir()
        (analysis/'analysis.complete').touch()
        merged=root/'mergedseq5';sintax=root/'sintax'
        hits=root/'best_hits.tsv'
        with hits.open('w') as handle:
            handle.write('merged_folder\tbarcode_folder\tseqID\tforward_patternName\treverse_patternName\n')
            for folder in ('barcode15','unclassified'):
                for raw,forward,reverse in [('r1','forward_tagged_primer_44','rev'),('r2','other','forward_tagged_primer_44'),('r3','bc1087',''),('r4','bc1087','')]:
                    handle.write(f'mergedseq5\t{folder}\t{raw}\t{forward}\t{reverse}\n')
        params=dict(merged_dir=str(merged),sintax_root=str(sintax),primer_hits=str(hits),cutoff=.8,analysed_folders=['barcode15','unclassified'])
        (analysis/'analysis_parameters.json').write_text(json.dumps(params))
        (analysis/'primer_name_mapping.json').write_text(json.dumps({'forward_tagged_primer_44':'bc1087','other':'bc1002','rev':'bc1100rc'}))
        (analysis/'taxon_plate_prevalence.tsv').write_text('rank\ttaxon_label\tcontamination\nspecies\tTipula maxima_123\tcontamination\n')
        rows=[]
        for folder in params['analysed_folders']:
            fastq(merged/folder/'raw.fastq.gz',[(raw,'AAAACCCCGGGG') for raw in ['r1','r2','r3','r4']])
            fastq(merged/folder/'filtered/asinput.fastq.gz',[
                ('split-by-adapter-left-r1','ACGT'),('fragment-x parent_read_id=r1','GGTT'),('r2','AAAA'),('r4','CCCC')])
            source=sintax/folder/f'{folder}_sintaxSP.tsv';source.parent.mkdir(parents=True)
            tax='g:Tipula_1(0.9),s:Tipula maxima_123(0.85)'
            source.write_text(f'split-by-adapter-left-r1\t{tax}\t+\t\nfragment-x parent_read_id=r1\t{tax}\t-\t\nr2\t{tax}\t+\t\nr4\t\t*\t\n')
            rows.append((folder,str(source.resolve()),source.stat().st_size,source.stat().st_mtime_ns))
        with (analysis/'input_manifest.tsv').open('w') as handle:
            handle.write('source_folder\tsintax_file\tsintax_bytes\tsintax_mtime_ns\n')
            for row in rows: handle.write('\t'.join(map(str,row))+'\n')
        return analysis

    def test_parent_trace_headers_and_unclassified_switch(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);analysis=self.fixture(root)
            out=root/'selected.fasta'
            self.assertEqual(export_primer(primer='bc1087',barcode='barcode15',analysis_dir=analysis,output=out),3)
            text=out.read_text()
            self.assertIn('>barcode15-bc1087-bc1100rc-Tipula_maxima_123-contamination-split-by-adapter-left-r1\nACGT\n',text)
            self.assertIn('>barcode15-bc1087-bc1100rc-Tipula_maxima_123-contamination-fragment-x\nGGTT\n',text)
            self.assertIn('>barcode15-bc1087-none-unassigned-nocontamination-r4\nCCCC\n',text)
            self.assertNotIn('unclassified-',text)
            self.assertNotIn('-r2\n',text)  # Reverse hits and nonwinning forward tags do not select reads.
            with out.with_suffix('.reads.tsv').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['raw_ont_read_id'] for r in rows if r['status']=='exported'],['r1','r1','r4'])
            self.assertEqual([r['raw_ont_read_id'] for r in rows if r['status']=='no_clean_read'],['r3'])
            self.assertTrue(all(r['forward_primers_detected']=='bc1087' for r in rows))
            self.assertTrue(all(r['raw_fastq_file'].endswith('/barcode15/raw.fastq.gz') for r in rows))
            taxout=root/'taxon.fasta'
            export_reads(taxa=['Tipula'],folders=['barcode15'],analysis_dir=analysis,output=taxout)
            self.assertTrue(set(taxout.read_text().splitlines()[::2]) - {'>barcode15-bc1002-bc1087-Tipula_maxima_123-contamination-r2'} <= set(text.splitlines()[::2]))
            both=root/'both.fasta'
            self.assertEqual(export_primer(primer='bc1087',barcode='barcode15',analysis_dir=analysis,output=both,include_unclassified=True),6)
            self.assertIn('>unclassified-',both.read_text())
            original=out.read_bytes()
            with self.assertRaises(FileExistsError): export_primer(primer='bc1087',barcode='barcode15',analysis_dir=analysis,output=out)
            self.assertEqual(out.read_bytes(),original)

    def test_low_confidence_species_uses_trusted_genus(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);analysis=self.fixture(root)
            manifest=analysis/'input_manifest.tsv'
            with manifest.open() as handle:
                reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
            for row in rows:
                source=Path(row['sintax_file']);source.write_text(source.read_text().replace('(0.85)','(0.79)'))
                row['sintax_bytes']=source.stat().st_size;row['sintax_mtime_ns']=source.stat().st_mtime_ns
            with manifest.open('w') as handle:
                writer=csv.DictWriter(handle,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
            output=root/'genus.fasta'
            export_primer(primer='bc1087',barcode='barcode15',analysis_dir=analysis,output=output)
            self.assertIn('barcode15-bc1087-bc1100rc-Tipula_1-nocontamination-',output.read_text())
            self.assertNotIn('Tipula_maxima',output.read_text())
            self.assertNotIn('taxapass',output.read_text())

    def test_queried_rank_confidence_independent_of_lower_ranks(self):
        for confidence in ('0.79', '0.80', '0.95', None):
            with self.subTest(confidence=confidence), tempfile.TemporaryDirectory() as d:
                root=Path(d);analysis=self.fixture(root)
                manifest=analysis/'input_manifest.tsv'
                with manifest.open() as handle:
                    reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
                for row in rows:
                    source=Path(row['sintax_file'])
                    support='' if confidence is None else f'({confidence})'
                    source.write_text(source.read_text().replace('g:Tipula_1(0.9)',f'f:Tipulidae_7147{support},g:Tipula_1(0.2)').replace('(0.85)','(0.1)'))
                    row['sintax_bytes']=source.stat().st_size;row['sintax_mtime_ns']=source.stat().st_mtime_ns
                with manifest.open('w') as handle:
                    writer=csv.DictWriter(handle,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
                for query in ('Tipulidae','7147','Tipulidae_7147'):
                    output=root/f'{query}.fasta'
                    if confidence in ('0.80','0.95'):
                        self.assertEqual(export_reads(taxa=[query],folders=['barcode15'],analysis_dir=analysis,output=output),3)
                        self.assertIn('-Tipulidae_7147-nocontamination-',output.read_text())
                    else:
                        with self.assertRaisesRegex(ValueError,'queried-rank confidence'):
                            export_reads(taxa=[query],folders=['barcode15'],analysis_dir=analysis,output=output)
                        self.assertFalse(output.exists())
                with self.assertRaisesRegex(ValueError,'queried-rank confidence'):
                    export_reads(taxa=['Tipula'],folders=['barcode15'],analysis_dir=analysis,output=root/'genus.fasta')

    def test_missing_raw_read_no_publication(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);analysis=self.fixture(root)
            fastq(root/'mergedseq5/barcode15/raw.fastq.gz',[('r1','AAAA'),('r2','AAAA'),('r4','AAAA')])
            with self.assertRaisesRegex(ValueError,'missing from raw'):
                export_primer(primer='bc1087',barcode='barcode15',analysis_dir=analysis,output=root/'fail.fasta')
            self.assertFalse(list(root.glob('fail.*')))

    def test_unanalysed_optional_folder_is_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);analysis=self.fixture(root)
            path=analysis/'analysis_parameters.json';params=json.loads(path.read_text());params['analysed_folders']=['barcode15'];path.write_text(json.dumps(params))
            with self.assertRaisesRegex(ValueError,'unclassified is not in this completed analysis'):
                export_primer(primer='bc1087',barcode='barcode15',analysis_dir=analysis,output=root/'fail.fasta',include_unclassified=True)
            self.assertFalse(list(root.glob('fail.*')))

if __name__=='__main__': unittest.main()
