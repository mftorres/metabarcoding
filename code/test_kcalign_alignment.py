"""Real installed kc-align integration, using temporary coding-sequence fixtures."""
import random
from pathlib import Path
import shutil
import tempfile
import unittest
from kcalign_alignment import align_kcalign, validate_reference
from kcalign_parallel import merge_results, screen_read
from unittest.mock import patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from postprocess_sintax import Runner, fasta

@unittest.skipUnless(shutil.which('kc-align'), 'kc-align unavailable')
class KcTests(unittest.TestCase):
    def test_real_mixed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);rng=random.Random(666)
            seq=''.join(rng.choice(['ATG','GCT','TTT','CCA','GGA']) for _ in range(80))
            reference=p/'ref.fasta';reference.write_text('>reference\n'+seq+'\n')
            inputs=p/'input.fasta'
            headers=['barcode01-bc1-bc2-Taxon-nocontamination-read1','REF1 description','BOLD_0001 description','REF2']
            inputs.write_text(''.join(f'>{h}\n{seq}\n' for h in headers))
            output=p/'combined_Test_nonumt_20261027.aln.fasta'
            try:
                align_kcalign(Runner(p),inputs,output,reference,threads=6)
            except Exception:
                print((p/'postprocessing.log').read_text())
                raise
            self.assertEqual(set(dict(fasta(output))),set(headers))
            self.assertTrue((p/'combined_Test_nonumt_20261027.pre_align.aa.fasta').is_file())
            self.assertFalse((p/'codon_align.fasta').exists())
    def test_rejected_read_keeps_id(self):
        ref={'guide':Seq('ATG'*80)}
        record=SeqRecord(Seq('ATG'*80),id='rejected_query')
        with patch('kcalign.find_homologs',return_value=1):
            result=screen_read((record,Seq('M'*80),ref,5,'none'))
        self.assertEqual(result[2],[])
        self.assertEqual(result[4],['rejected_query'])
        accepted=([Seq('M')],['description'],['accepted'],{'accepted':Seq('ATG')},[])
        combined=merge_results([result,accepted],ref)
        self.assertEqual(combined[2],['accepted'])
        self.assertEqual(combined[4],['rejected_query'])
        self.assertNotIn('rejected_query',combined[3])

    def test_invalid_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'reference.fasta'
            for text in ('>x\nAT\n','>x\nATN\n','>x\nATG\n>y\nATG\n','>x\nTAAGCT\n'):
                p.write_text(text)
                with self.assertRaises(ValueError):validate_reference(p)

if __name__=='__main__':unittest.main()
