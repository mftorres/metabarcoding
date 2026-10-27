"""Isolated tests for confidence, plate prevalence and clean-read export."""
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
import pandas as pd
from sintax_parser import parse_output
from sintax_analysis import run_analysis, prevalence_table
from export_sintax_fasta import export_reads


class ReadAnalysisTests(unittest.TestCase):
    def test_rank_assignment_at_boundary(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'sample_sintaxSP.tsv'
            p.write_text('r\tc:Insecta_50557(1.00),g:Fly_1(0.80),s:Fly one_2(0.79)\t+\tc:Insecta_50557\n')
            row=next(parse_output(p,cutoff=.8))
            self.assertEqual(row['assigned_rank'],'genus')
            self.assertEqual(row['assigned_confidence'],.8)
            self.assertEqual(row['taxonomy_status'],'taxapass')
            self.assertEqual(row['taxlabel'],'c:Insecta_50557')
            self.assertEqual(row['species_bootstrap'],.79)

    def test_distinct_plates_and_eighty_percent(self):
        p=prevalence_table({('species','Fly_2'):{'a','b','c','d'},('genus','Fly_1'):{'a','b','c'}},set('abcde'),.8)
        self.assertEqual(p.set_index('rank').loc['species','contamination'],'contamination')
        self.assertEqual(p.set_index('rank').loc['genus','contamination'],'nocontamination')
        self.assertEqual(p.minimum_plate_count.tolist(),[4,4])

    def test_end_to_end_analysis_and_clean_export(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);data=root/'data';data.mkdir()
            (data/'saproxdna_primertags_seq5_2026_1.fasta').write_text('>f1\nACGT\n>r1\nTGCA\n')
            (data/'PCRsaproxdna_primertags_seq5_2026_1.fasta').write_text('>bc1001\nACGT\n>bc1100rc\nTGCA\n')
            merged=data/'mergedseq5';sintax=data/'sintax';out=data/'analysis'
            folders=[f'barcode{i:02d}' for i in range(1,7)]
            primer_rows=[]
            taxonomy='c:Insecta_50557(1.00),o:Diptera_10(1.00),f:Family_11(1.00),g:Fly_1(0.95),s:Fly one_2(0.90)'
            for i,folder in enumerate(folders[:5]):
                sample=merged/folder/'filtered';sample.mkdir(parents=True)
                with gzip.open(sample/'asinput.fastq.gz','wt') as handle:
                    handle.write('@split-by-adapter-left-r1\nACGTAC\n+\nIIIIII\n@r2\nTTGGCC\n+\nIIIIII\n')
                p=sintax/folder/f'{folder}_sintaxSP.tsv';p.parent.mkdir(parents=True)
                text=taxonomy if i<4 else taxonomy.replace('(0.90)','(0.79)')
                p.write_text(f'split-by-adapter-left-r1\t{text}\t+\tc:Insecta_50557\nr2\t{text.replace("(0.90)","(0.80)")}\t-\t\n')
                primer_rows.append({'merged_folder':'mergedseq5','barcode_folder':folder,'seqID':'r1','forward_patternName':'f1','reverse_patternName':'r1'})
                primer_rows.append({'merged_folder':'mergedseq5','barcode_folder':folder,'seqID':'r2','forward_patternName':'','reverse_patternName':'r1'})
            primer=root/'primers.tsv';pd.DataFrame(primer_rows).to_csv(primer,sep='\t',index=False)
            plates=root/'plates.tsv';pd.DataFrame({'sequencing_barcode':folders,'PCR_plate':[f'p{i}' for i in range(6)]}).to_csv(plates,sep='\t',index=False)
            result=run_analysis(project_root=root,merged_dir=merged,sintax_root=sintax,output_dir=out,
                primer_hits=primer,plate_map=plates,folders=folders,allow_missing=True)
            self.assertEqual(result['analysed_plate_count'],5)
            self.assertEqual(result['missing_folders'],['barcode06'])
            summary=pd.read_csv(out/'all_species_summary.tsv',sep='\t')
            self.assertEqual(summary.read_count.tolist(),[2,2,2,2,0])
            self.assertAlmostEqual(summary.median_confidence.iloc[0],.85)
            self.assertEqual(summary.forward_primers_detected.iloc[0],'bc1001')
            self.assertEqual(summary.reverse_primers_detected.iloc[0],'bc1100rc')
            self.assertNotIn('forward_primer_names_original', summary.columns)
            self.assertEqual(summary['class'].iloc[0], 'Insecta')
            prevalence=pd.read_csv(out/'taxon_plate_prevalence.tsv',sep='\t')
            self.assertTrue(prevalence.reverse_primers_detected.eq('bc1100rc').all())
            pairs=pd.read_csv(out/'barcode01/barcode01_primer_taxonomy.tsv',sep='\t',keep_default_na=False)
            self.assertTrue((pairs.forward_PCR_primer.eq('') & pairs.reverse_PCR_primer.eq('bc1100rc')).any())
            self.assertEqual(summary.contamination.tolist(),['contamination']*5)
            self.assertTrue(pd.isna(summary.median_confidence.iloc[-1]))
            reads=pd.read_csv(out/'all_reads.tsv',sep='\t')
            self.assertEqual(len(reads),10)
            self.assertNotIn('forward_patternName', reads.columns)
            self.assertNotIn('reverse_patternName', reads.columns)
            self.assertTrue(reads.reverse_PCR_primer.eq('bc1100rc').all())
            self.assertEqual(reads.source_folder.nunique(),5)
            self.assertEqual(reads.raw_ont_read_id.iloc[0],'r1')
            fasta=root/'selected.fasta'
            count=export_reads(taxa=['insecta'],folders=['barcode01','barcode05'],analysis_dir=out,output=fasta)
            self.assertEqual(count,4)
            text=fasta.read_text()
            self.assertIn('>barcode01-bc1001-bc1100rc-Fly_one_2-contamination-split-by-adapter-left-r1\nACGTAC',text)
            self.assertIn('>barcode05-none-bc1100rc-Fly_1-contamination-r2',text)
            with self.assertRaises(FileExistsError):
                export_reads(taxa=['Insecta'],folders=['barcode01'],analysis_dir=out,output=fasta)
            with self.assertRaises(ValueError):
                export_reads(taxa=['no-such-taxon'],folders=['barcode01'],analysis_dir=out,output=root/'missing.fasta')
            self.assertFalse((root/'missing.fasta').exists())
            with self.assertRaises(ValueError):
                export_reads(taxa=['Insecta'],folders=['barcode06'],analysis_dir=out,output=root/'missing.fasta')

if __name__=='__main__': unittest.main()
