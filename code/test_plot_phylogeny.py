"""Real Baltic rendering and support-preservation tests using temporary fixtures."""
import csv
import itertools
from pathlib import Path
import tempfile
import unittest
import baltic as bt
from plot_phylogeny import render, midpoint_with_support, split_keys


def distances(tree):
    def path(leaf):
        result={leaf:0.0}
        distance=0
        while leaf is not tree.root:
            distance+=leaf.length
            leaf=leaf.parent
            result[leaf]=distance
        return result
    leaves={b.name:b for b in tree.get_external()}
    result={}
    for a,b in itertools.combinations(sorted(leaves),2):
        pa,pb=path(leaves[a]),path(leaves[b])
        result[a,b]=min(pa[k]+pb[k] for k in pa.keys() & pb.keys())
    return result


class PlotTests(unittest.TestCase):
    def test_midpoint_support_and_distances(self):
        tree=bt.make_tree('((A:0.1,B:0.2)91/100:0.3,(C:0.15,D:0.2)88:0.1,E:0.7);',treeType='divergence')
        before=distances(tree)
        expected={key:b.traits['label'] for b,key in split_keys(tree).items() if 'label' in b.traits}
        midpoint_with_support(tree)
        after=distances(tree)
        self.assertEqual(expected,{key:b.traits['label'] for b,key in split_keys(tree).items() if 'label' in b.traits})
        for key in before: self.assertAlmostEqual(before[key],after[key])
        farthest=max(before,key=before.get)
        def depth(name):
            branch=next(b for b in tree.get_external() if b.name==name)
            value=0
            while branch is not tree.root:
                value+=branch.length;branch=branch.parent
            return value
        self.assertAlmostEqual(depth(farthest[0]),before[farthest]/2)
        self.assertAlmostEqual(depth(farthest[1]),before[farthest]/2)

    def test_figure_and_labels(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            folder=p/'analysis/barcode01';folder.mkdir(parents=True)
            local='barcode01-bc1087-bc1100rc-Tipula_123-nocontamination-read-uuid'
            with (folder/'barcode01_reads.tsv').open('w') as f:
                f.write('readid\taccepted_species\tassigned_name\nread-uuid\tTipula_maxima_123\tTipula maxima\n')
            (p/'ref.xml').write_text('<GBSet>'+''.join(f'<GBSeq><GBSeq_accession-version>{a}</GBSeq_accession-version><GBSeq_organism>Tipula test{i}</GBSeq_organism></GBSeq>' for i,a in enumerate(['REF1.1','REF2.1','REF3.1']))+'</GBSet>')
            (p/'tree').write_text(f'(({local}:0.1,REF1.1:0.2)95/99:0.3,REF2.1:0.2,REF3.1:0.7);')
            render(p/'tree',p/'analysis',p/'ref.xml',p/'figure')
            self.assertGreater((p/'figure.pdf').stat().st_size,1000)
            self.assertIn('subs/site',(p/'figure.svg').read_text())
            with (p/'figure.tips.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertTrue(any(r['display_label']=='Tipula maxima | read-uuid | barcode01' for r in rows))
            self.assertTrue(any(r['display_label']=='Tipula test0 | REF1.1' for r in rows))
            self.assertIn('95/99',(p/'figure.svg').read_text())
            (p/'run.treefile').write_text((p/'tree').read_text().replace('95/99',''))
            (p/'run.contree').write_text((p/'tree').read_text())
            render(p/'run.treefile',p/'analysis',p/'ref.xml',p/'consensus')
            self.assertIn('run.contree',(p/'consensus.plot.json').read_text())
            with self.assertRaises(FileExistsError): render(p/'tree',p/'analysis',p/'ref.xml',p/'figure')

if __name__=='__main__':unittest.main()
