"""Real VSEARCH tests for local grouping and global reference clustering."""
import csv
from pathlib import Path
import random
import shutil
import tempfile
import unittest
from postprocess_sintax import Runner, cluster_local_groups, fasta, vsearch_cluster

@unittest.skipUnless(shutil.which('vsearch'), 'VSEARCH unavailable')
class VsearchGroupingTests(unittest.TestCase):
    def test_group_isolation_and_reference_consensus(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);runner=Runner(root)
            rng=random.Random(666);sequence=''.join(rng.choice('ACGT') for _ in range(650))
            ids=['barcode01-bc1001-bc1100rc-Fly_1-nocontamination-read1',
                 'barcode01-bc1001-bc1100rc-Fly_1-nocontamination-read2',
                 'barcode01-bc1002-bc1100rc-Fly_1-nocontamination-read3',
                 'barcode01-bc1001-bc1101rc-Fly_1-nocontamination-read4',
                 'barcode02-bc1001-bc1100rc-Fly_1-nocontamination-read5',
                 'barcode01-bc1001-bc1100rc-Other_2-nocontamination-read6']
            source=root/'local.fasta';source.write_text(''.join(f'>{name}\n{sequence}\n' for name in ids))
            original=source.read_bytes()
            cluster_local_groups(runner,source,root/'centroids.fasta',root/'groups',root/'groups.tsv',root/'memberships.tsv',1)
            records=list(fasta(root/'centroids.fasta'))
            self.assertEqual(len(records),5)
            self.assertEqual({tuple(h.split('-')[:4]) for h,s in records},{tuple(h.split('-')[:4]) for h in ids})
            self.assertTrue(all(h in ids and s==sequence for h,s in records))
            self.assertEqual(source.read_bytes(),original)
            with (root/'groups.tsv').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(sorted(int(row['input_sequences']) for row in rows),[1,1,1,1,2])
            with (root/'memberships.tsv').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len([r for r in rows if r['record_type'] in ('S','H')]),6)
            refs=root/'refs.fasta';refs.write_text(f'>REF1.1 species one\n{sequence}\n>REF2.1 species two\n{sequence}\n')
            before,after=vsearch_cluster(runner,refs,root/'refs.centroids.fasta',root/'refs.uc',root/'refs.log',1)
            self.assertEqual((before,after),(2,1))
            self.assertIn(list(fasta(root/'refs.centroids.fasta'))[0][0],['REF1.1 species one','REF2.1 species two'])


    def test_consensus_differs_from_seed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); runner=Runner(root)
            rng=random.Random(42)
            majority=''.join(rng.choice('ACGT') for _ in range(650))
            seed=majority[:300]+('A' if majority[300]!='A' else 'C')+majority[301:]
            source=root/'input.fasta'
            source.write_text(f'>aaa_seed description\n{seed}\n>member2\n{majority}\n>member3\n{majority}\n')
            self.assertEqual(vsearch_cluster(runner,source,root/'consensus.fasta',root/'clusters.uc',root/'stdout',1),(3,1))
            records=list(fasta(root/'consensus.fasta'))
            self.assertEqual(records,[('aaa_seed description',majority)])
            self.assertNotEqual(records[0][1],seed)

if __name__=='__main__':unittest.main()
