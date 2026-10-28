"""Validate installed ClipKIT smart-gap behavior on an isolated nucleotide alignment."""
import subprocess
import tempfile
from pathlib import Path
import shutil
import unittest
from postprocess_sintax import fasta

@unittest.skipUnless(shutil.which('clipkit'), 'ClipKIT unavailable')
class ClipkitTests(unittest.TestCase):
    def test_smart_gap_keeps_headers_and_removes_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.fasta';output=root/'pruned.fasta'
            records=[('read1 description','ACGT--ACGT'),('read2','ACGT--ACGT'),('read3','ACGT--ACGT'),('read4','ACGT--ACGT')]
            source.write_text(''.join(f'>{h}\n{s}\n' for h,s in records))
            before=source.read_bytes()
            subprocess.run(['clipkit',str(source),'-m','smart-gap','-s','nt','-if','fasta','-of','fasta','-o',str(output)],check=True,capture_output=True)
            self.assertEqual(list(fasta(output)),[(h,'ACGTACGT') for h,_ in records])
            self.assertEqual(source.read_bytes(),before)

if __name__=='__main__':unittest.main()
