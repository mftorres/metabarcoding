"""Exercise installed MACSE on temporary coding sequences and validate resume hashes."""
import argparse
import json
from pathlib import Path
import random
import shutil
import tempfile
import unittest
from macse_alignment import align_macse
from postprocess_sintax import Runner, fasta, file_digest, latest_sequences

@unittest.skipUnless(shutil.which('macse'), 'MACSE unavailable')
class MacseTests(unittest.TestCase):
    def test_real_alignment_trim_and_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            retained=root/'combined_Test_nonumt_20261025.fasta'
            alignment=retained.with_suffix('.aln.fasta')
            rng=random.Random(4)
            seq=''.join(rng.choice(['ATG','GCT','TTT','CCA','GGA']) for _ in range(80))
            headers=['barcode01-bc1-bc2-Test-nocontamination-read1', 'REF1 description', 'BOLD_000000001 original header', 'REF2']
            retained.write_text(''.join(f'>{h}\n{seq}\n' for h in headers))
            align_macse(Runner(root), retained, alignment, memory="512m")
            self.assertEqual(set(dict(fasta(alignment))),set(headers))
            self.assertTrue(all('!' not in seq for _,seq in fasta(alignment)))
            self.assertTrue((root/'combined_Test_nonumt_20261025_trimmed.csv').is_file())
            marker=root/'alignment.complete.json'
            marker.write_text(json.dumps(dict(algorithm='MACSE',alignment_sha256=file_digest(alignment),retained_sha256=file_digest(retained))))
            args=argparse.Namespace(date='20261025',macse=True,skip_translation=False,min_aa=200,aa_length_mode='stop-free')
            self.assertTrue(latest_sequences(root,args,{})[2])
            alignment.write_text(alignment.read_text().replace('ATG','GCT',1))
            self.assertFalse(latest_sequences(root,args,{})[2])

if __name__=='__main__': unittest.main()
