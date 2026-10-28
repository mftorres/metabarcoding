"""Tests use temporary fixtures; no existing analysis files or network."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import postprocess_sintax as m

class Tests(unittest.TestCase):
    def test_reference_policy(self):
        self.assertIn(' NOT (', m.reference_query('7147'))
        self.assertNotIn(' NOT (', m.reference_query('7147', True))
        self.assertIn('COI[Gene Name]', m.reference_query('7147', True))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'download').write_text('>COI.1 mitochondrion cytochrome oxidase I gene, complete cds\nATG\n'
                '>GENOME.1 mitochondrion, complete genome\nATG\n'
                '>PART.1 partial mitochondrial genome\nATG\n'
                '>MITO.1 mitogenome\nATG\n'
                '>SEQ.1 mitochondrion, complete sequence\nATG\n')
            self.assertEqual(m.filter_references(p/'download',p/'out',p/'audit',p/'ids'),(1,4))
            self.assertEqual((p/'ids').read_text(),'COI.1\n')
            self.assertEqual(m.filter_references(p/'download',p/'all',p/'auditall',p/'idsall',True),(5,0))
            (p/'empty').write_text('>GENOME.1 complete genome\nATG\n')
            with self.assertRaisesRegex(ValueError,'No references remain'):
                m.filter_references(p/'empty',p/'none',p/'auditnone',p/'idsnone')

    def test_filter(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'nt').write_text('>pass original header\nATG\n>fail\nATG\n')
            with (p/'aa').open('w') as f:
                for name in ['pass', 'fail']:
                    for frame in [1,2,3,-1,-2,-3]:
                        m.write_record(f, f'{name}_frame={frame}', 'M'*200 if name=='pass' and frame==-2 else 'M'*199+'*M')
            self.assertEqual(m.filter_translations(p/'nt',p/'aa',p/'out',p/'report',200,'stop-free'),1)
            self.assertEqual(list(m.fasta(p/'out')),[('pass original header','ATG')])
            self.assertEqual(m.filter_translations(p/'nt',p/'aa',p/'total',p/'report2',200,'total'),2)

    @unittest.skipUnless(shutil.which('seqkit'), 'seqkit unavailable')
    def test_real_seqkit(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'nt').write_text('>read description\n'+'ATG'*200+'\n')
            with (p/'aa').open('w') as f:
                subprocess.run(['seqkit','translate','-f','6','-F','-T','5',str(p/'nt')],stdout=f,check=True)
            self.assertEqual(len(list(m.fasta(p/'aa'))),6)
            self.assertEqual(m.filter_translations(p/'nt',p/'aa',p/'out',p/'report',200,'stop-free'),1)

    @unittest.skipUnless(shutil.which('seqkit'), 'seqkit unavailable')
    def test_reference_length_boundaries(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'input').write_text(''.join(f'>n{n}\n'+ 'A'*n+'\n' for n in (499,500,2000,2001)))
            with (p/'output').open('w') as f:
                subprocess.run(['seqkit','seq','--min-len','500','--max-len','2000',str(p/'input')],stdout=f,check=True)
            self.assertEqual([len(seq) for h,seq in m.fasta(p/'output')],[500,2000])

    def test_workflow(self):
        self.workflow("download")

    def test_bold_workflow(self):
        self.workflow("bold")

    def test_skip_translation_workflow(self):
        self.workflow('bold', skip=True)

    def test_kcalign_workflow(self):
        self.workflow('bold', kcalign=True)

    def test_kcalign_without_translation_workflow(self):
        self.workflow('bold', skip=True, kcalign=True)

    def workflow(self, source, skip=False, kcalign=False):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'analysis.complete').touch()
            args=argparse.Namespace(taxon='Tipulidae',folders=['barcode01','barcode02'],analysis_dir=p,output_dir=p/'run',date='20261021',threads=6,translation_table=5,min_aa=200,aa_length_mode='stop-free',bootstrap_method='standard',include_mitogenomes=False)
            args.kcalign=kcalign
            if kcalign:
                args.kcalign_reference=p/'trusted.fasta'
                args.kcalign_reference.write_text('>trusted\n'+'ATG'*200+'\n')
            args.skip_translation=skip
            args.reference_source=source
            if source == 'bold':
                import gzip
                args.bold_db=p/'bold.fasta.gz'
                with gzip.open(args.bold_db,'wt') as handle:
                    handle.write('>PROCESS-1|COI-5P|UK|Animalia,Tipulidae,Tipula maxima,None\n'+'ATG'*200+'\n')
            commands=[]
            fail_alignment=[True]
            def fake(runner, command, output, input_path=None, cwd=None):
                c=list(map(str,command)); commands.append(c); tool=Path(c[0]).name
                if tool=='esearch': Path(output).write_text('<ENTREZ_DIRECT/>')
                elif tool=='efetch':
                    fmt=c[c.index('-format')+1]
                    Path(output).write_text({'xml':'<TaxaSet><Taxon><TaxId>7147</TaxId><ScientificName>Tipulidae</ScientificName></Taxon></TaxaSet>','gb':'<GBSet/>','acc':'REF.1\n','fasta':'>REF.1 description\n'+'ATG'*200+'\n'}[fmt])
                elif len(c)>1 and c[1].endswith('export_sintax_taxa.py'):
                    Path(c[c.index('--output')+1]).write_text(''.join(f'>{folder}-bc1087-none-Tipula-nocontamination-read{i}\n'+'ATG'*200+'\n' for i,folder in enumerate(['barcode01','barcode01','barcode02'])))
                    Path(output).write_text('Exported 3 reads\n')
                elif len(c)>1 and c[1].endswith('plot_phylogeny.py'):
                    Path(output).write_text('Mock figure invocation\n')
                    for suffix in ('.pdf','.svg','.midpoint.nwk','.tips.tsv','.plot.json'):
                        Path(c[c.index('--output-prefix')+1]+suffix).write_text('mock figure output\n')
                elif tool=='vsearch':
                    Path(c[c.index('--consout')+1]).write_text(''.join(f'>centroid={h};seqs=1\n{seq}\n' for h,seq in m.fasta(c[c.index('--cluster_fast')+1]))); Path(output).touch()
                    Path(c[c.index('--uc')+1]).write_text('S\t0\t600\t*\t*\t*\t*\t*\tread\t*\n')
                elif tool=='seqkit' and c[1]=='seq':
                    with Path(output).open('w') as f:
                        for h,seq in m.fasta(c[-1]):
                            if 500 <= len(seq) <= 2000: m.write_record(f,h,seq)
                elif tool=='seqkit':
                    with Path(output).open('w') as f:
                        for h,s in m.fasta(c[-1]):
                            for frame in [1,2,3,-1,-2,-3]: m.write_record(f,h.split()[0]+f'_frame={frame}','M'*200)
                elif len(c)>1 and c[1].endswith('kcalign_cli.py'):
                    Path(output).write_text('kc-align log\n')
                    if fail_alignment[0]:
                        fail_alignment[0]=False
                        raise subprocess.CalledProcessError(1,c)
                    data=Path(c[c.index('--sequences')+1]).read_text()+Path(c[c.index('--reference')+1]).read_text()
                    (Path(cwd)/'codon_align.fasta').write_text(data)
                    (Path(cwd)/'pre_align.fasta').write_text(data)
                elif tool=='mafft':
                    if fail_alignment[0]:
                        fail_alignment[0]=False
                        Path(output).write_text('partial alignment\n')
                        raise subprocess.CalledProcessError(1,c)
                    with Path(output).open('w') as f:
                        for i,(h,s) in enumerate(m.fasta(c[-1])): m.write_record(f,('_R_' if i==0 else '')+h,s)
                elif tool=='clipkit':
                    shutil.copyfile(c[1],c[c.index('-o')+1]);Path(output).write_text('trimmed fixture\n')
                elif tool=='iqtree':
                    Path(c[c.index('--prefix')+1]+'.treefile').write_text('(a,b,c,d);\n'); Path(output).touch()
                else: raise AssertionError(c)
            with patch.object(m.shutil,'which',side_effect=lambda n:'/mock/'+n),patch.object(m.Runner,'run',fake):
                with self.assertRaises(subprocess.CalledProcessError): m.run_pipeline(args)
                retained=next((p/'run').glob('*selected*.fasta' if skip else '*nonumt*.fasta'))
                retained_digest=m.file_digest(retained)
                self.assertFalse(list((p/'run').glob('*.aln.fasta*')))
                # Existing batched metadata must never trigger another efetch on resume.
                (p/'run/reference_records.gb.xml').write_text('<?xml version="1.0"?><GBSet/>\n<?xml version="1.0"?><GBSet/>')
                args.resume=True
                args.date=None
                prior=len(commands)
                m.run_pipeline(args)
                resumed=commands[prior:]
                self.assertEqual([Path(c[0]).name for c in resumed[:2]],[Path(m.sys.executable).name,'clipkit'] if kcalign else ['mafft','clipkit'])
                self.assertEqual(len(resumed),4)
                self.assertEqual(m.file_digest(retained),retained_digest)
                self.assertTrue((p/'run/alignment.complete.json').exists())
                self.assertFalse(any(c[0]=='efetch' for c in resumed))
                self.assertFalse((p/'run/postprocessing.failed').exists())
                prior=len(commands)
                m.run_pipeline(args)
                self.assertEqual(len(commands),prior)
                # Explicit restart reuses retained sequences and rebuilds only alignment/tree/figure.
                args.restart_alignment=True
                prior=len(commands)
                m.run_pipeline(args)
                self.assertEqual(len(commands)-prior,4)
                self.assertEqual(m.file_digest(retained),retained_digest)
                self.assertFalse(any(c[0] in ('vsearch','seqkit','efetch','esearch') for c in commands[prior:]))
                args.restart_alignment=False
                prior=len(commands)
                m.run_pipeline(args)
                self.assertEqual(len(commands),prior)
                # Completed alignment survives a missing pipeline checkpoint: start at IQ-TREE.
                (p/'run/checkpoints.json').unlink()
                prior=len(commands)
                m.run_pipeline(args)
                self.assertEqual(Path(commands[prior][0]).name,'clipkit')
                self.assertEqual(len(commands[prior:]),3)
                self.assertEqual(m.file_digest(retained),retained_digest)
                if kcalign:
                    command=next(c for c in commands if len(c)>1 and c[1].endswith('kcalign_cli.py'))
                    self.assertEqual(command[command.index('--threads')+1],'6')
                    self.assertIn('--keep',command)
                    self.assertFalse(any(c[0]=='mafft' for c in commands))
                else:
                    mafft=next(c for c in commands if c[0]=='mafft')
                    self.assertEqual(mafft[mafft.index('--retree')+1],'2')
                    self.assertEqual(mafft[mafft.index('--maxiterate')+1],'0')
                    self.assertNotIn('--localpair',mafft)
                import json
                parameter_path=p/'run/parameters.json'
                saved_parameters=parameter_path.read_text()
                old_parameters=json.loads(saved_parameters);old_parameters.pop('header_schema')
                parameter_path.write_text(json.dumps(old_parameters))
                with self.assertRaisesRegex(ValueError,'old FASTA headers'): m.run_pipeline(args)
                parameter_path.write_text(saved_parameters)
                for option in ('skip_translation',):
                    original=getattr(args,option);setattr(args,option,not original)
                    with self.assertRaisesRegex(ValueError,'changed '+option): m.run_pipeline(args)
                    setattr(args,option,original)
                if kcalign:
                    original=args.kcalign_reference.read_text()
                    args.kcalign_reference.write_text(original.replace('ATG','GCT'))
                    with self.assertRaisesRegex(ValueError,'changed trusted alignment reference'): m.run_pipeline(args)
                    args.kcalign_reference.write_text(original)
                args.taxon='Other'
                with self.assertRaisesRegex(ValueError,'changed taxon'): m.run_pipeline(args)
                args.taxon='Tipulidae'
                args.resume=False
                with self.assertRaises(FileExistsError): m.run_pipeline(args)
            if source == 'bold':
                self.assertFalse(any(Path(c[0]).name in ('efetch','esearch') for c in commands))
                self.assertTrue((p/'run/bold_reference_metadata.tsv').exists())
                args.resume=True
                args.reference_source='download'
                with patch.object(m.shutil,'which',return_value='/mock/tool'):
                    with self.assertRaisesRegex(ValueError,'Cannot change --download/--bold'): m.run_pipeline(args)
                args.reference_source='bold'
                args.bold_db.touch()
                with patch.object(m.shutil,'which',return_value='/mock/tool'):
                    with self.assertRaisesRegex(ValueError,'BOLD database files changed'): m.run_pipeline(args)
            self.assertTrue((p/'run/postprocessing.complete').exists())
            self.assertIn('barcode01\t2',(p/'run/extraction_counts.tsv').read_text())
            self.assertIn('barcode02\t1',(p/'run/extraction_counts.tsv').read_text())
            if skip:
                self.assertFalse(any(Path(c[0]).name=='seqkit' and c[1]=='translate' for c in commands))
                self.assertFalse((p/'run/translation_lengths.tsv').exists())
            clusters=[c for c in commands if c[0]=='vsearch']
            self.assertEqual(len(clusters),3)
            self.assertTrue(all(c[c.index('--id')+1]=='0.98' for c in clusters))
            inputs=[c[c.index('--cluster_fast')+1] for c in clusters]
            self.assertTrue(Path(inputs[0]).name.endswith('.input.fasta'))
            self.assertTrue(inputs[-1].endswith('.length500-2000.fasta'))
            self.assertEqual(len(list(m.fasta(next((p/'run').glob('combined_*.fasta'))))),4)
            c=next(c for c in commands if c[0]=='iqtree')
            self.assertTrue(c[c.index('-s')+1].endswith('.clipkit.aln.fasta'))
            self.assertEqual(c[c.index('-b')+1],'1000')
            self.assertEqual(c[c.index('-m')+1],'MFP')
            self.assertFalse(any(h.startswith('_R_') for h,s in m.fasta(next((p/'run').glob('*.aln.fasta')))))

    def test_checkpoint_tampering_invalidates_downstream(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            runner=m.Runner(p)
            checkpoints=m.Checkpoints(runner,False)
            for name in ('first','second'):
                with checkpoints.stage(name,[p/name]) as execute:
                    self.assertTrue(execute)
                    (p/name).write_text('original')
            (p/'first').write_text('modified')
            checkpoints=m.Checkpoints(runner,True)
            with checkpoints.stage('first',[p/'first']) as execute:
                self.assertTrue(execute)
                (p/'first').write_text('rebuilt')
            # Simulate interruption before rebuilding the downstream stage.
            checkpoints=m.Checkpoints(runner,True)
            with checkpoints.stage('first',[p/'first']) as execute: self.assertFalse(execute)
            with checkpoints.stage('second',[p/'second']) as execute:
                self.assertTrue(execute)
                (p/'second').write_text('rebuilt')
            self.assertTrue(list((p/'resume_archive').glob('*first/first')))

    @unittest.skipUnless(shutil.which('mafft'), 'mafft unavailable')
    def test_real_fft_ns2(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'input.fasta').write_text('>a\nACGTACGTACGT\n>b\nACGTACGTACGA\n>c\nACGTACGTACGG\n>d\nACGTACGTACGC\n')
            with (p/'aligned.fasta').open('w') as handle:
                subprocess.run(['mafft','--retree','2','--maxiterate','0','--adjustdirectionaccurately',
                                '--thread','1',str(p/'input.fasta')],stdout=handle,stderr=subprocess.PIPE,check=True)
            records=list(m.fasta(p/'aligned.fasta'))
            self.assertEqual(len(records),4)
            self.assertEqual(len({len(seq) for h,seq in records}),1)

    def test_missing_dependency(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'out'
            with patch.object(m.shutil,'which',return_value=None):
                with self.assertRaisesRegex(ValueError,'Missing executables'): m.run_pipeline(argparse.Namespace(output_dir=p))
            self.assertFalse(p.exists())

if __name__=='__main__': unittest.main()
