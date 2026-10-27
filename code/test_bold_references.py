"""BOLD extraction and actual Baltic rendering with isolated temporary files."""
import csv
import gzip
from pathlib import Path
import tempfile
import unittest
from bold_references import database_files, extract_bold
from plot_phylogeny import render


class BoldTests(unittest.TestCase):
    def test_selection_provenance_and_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            with gzip.open(p/'db.fasta.gz', 'wt') as handle:
                handle.write('>OTHER|COI-5P|UK|Animalia,Culicidae\nACGT\n')
                for i in range(4):
                    handle.write(f'>PROC-{i}|COI-5P|UK|Animalia,Tipulidae,Tipula maxima,None\nACG\nTAA\n')
            self.assertEqual(extract_bold(database_files(p), 'tipULIDAE', p/'out', p/'metadata.tsv', p/'ids'),4)
            with (p/'metadata.tsv').open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(rows[0]['taxon'], 'Tipula maxima')
            self.assertEqual(rows[0]['process_id'], 'PROC-0')
            self.assertIn('ACGTAA', (p/'out').read_text())
            self.assertNotIn('OTHER', (p/'out').read_text())
            self.assertEqual(extract_bold([p/'db.fasta.gz'], 'Tipula_maxima', p/'species', p/'species.tsv', p/'species.ids'),4)
            ids = (p/'ids').read_text().splitlines()
            (p/'tree').write_text(f'(({ids[0]}:0.1,{ids[1]}:0.2)95:0.3,{ids[2]}:0.2,{ids[3]}:0.7);')
            render(p/'tree', p, p/'metadata.tsv', p/'figure')
            with (p/'figure.tips.tsv').open() as handle:
                tips = list(csv.DictReader(handle, delimiter='\t'))
            self.assertTrue(all(t['source']=='bold' for t in tips))
            self.assertTrue(any(t['display_label']=='Tipula maxima | PROC-0' for t in tips))
            self.assertGreater((p/'figure.pdf').stat().st_size,1000)

    def test_reject_numeric_and_no_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            with gzip.open(p/'db.fa.gz', 'wt') as handle:
                handle.write('>PROC|COI|UK|Animalia,Tipulidae\nATG\n')
            with self.assertRaisesRegex(ValueError, 'numeric taxid'):
                extract_bold([p/'db.fa.gz'], '7147', p/'out', p/'meta', p/'ids')
            with self.assertRaisesRegex(ValueError, 'No BOLD FASTA headers'):
                extract_bold([p/'db.fa.gz'], 'Coleoptera', p/'out', p/'meta', p/'ids')

if __name__ == '__main__':
    unittest.main()
