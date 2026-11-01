"""Reference taxon suffixes must survive filtering and resolve during plotting."""
import csv
from pathlib import Path
import tempfile
import unittest
from postprocess_sintax import filter_references, fasta
from plot_phylogeny import render
from reference_labels import labelled_reference_id, reference_id


class ReferenceLabelTests(unittest.TestCase):
    def test_id_is_safe_and_idempotent(self):
        label = labelled_reference_id('REF.1', 'Tipula (Tipula) maxima')
        self.assertEqual(label, 'REF.1__taxon__Tipula_Tipula_maxima')
        self.assertEqual(reference_id(label), 'REF.1')
        self.assertEqual(labelled_reference_id(label, 'Tipula maxima'), 'REF.1__taxon__Tipula_maxima')

    def test_selection_and_plotting(self):
        for source in ('bold', 'ncbi'):
            with self.subTest(source=source), tempfile.TemporaryDirectory() as d:
                p = Path(d)
                taxa = {f'REF{i}.1': f'Tipula species {i}' for i in range(4)}
                (p/'download').write_text(''.join(f'>{key} original description\nACGT\n' for key in taxa))
                filter_references(p/'download', p/'selected', p/'audit', p/'ids', reference_taxa=taxa)
                selected = list(fasta(p/'selected'))
                self.assertEqual([seq for _, seq in selected], ['ACGT'] * 4)
                self.assertEqual((p/'ids').read_text().splitlines(), list(taxa))
                self.assertTrue(all(not any(c.isspace() for c in header) for header, _ in selected))
                a,b,c,d = [header for header, _ in selected]
                (p/'tree').write_text(f'(({a}:0.1,{b}:0.2)95:0.3,{c}:0.2,{d}:0.7);')
                if source == 'bold':
                    metadata = p/'metadata.tsv'
                    with metadata.open('w') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(['sequence_id', 'taxon', 'process_id'])
                        writer.writerows((key, value, f'PROC-{i}') for i,(key,value) in enumerate(taxa.items()))
                else:
                    metadata = p/'metadata.xml'
                    metadata.write_text('<GBSet>' + ''.join(f'<GBSeq><GBSeq_accession-version>{key}</GBSeq_accession-version><GBSeq_organism>{value}</GBSeq_organism></GBSeq>' for key,value in taxa.items()) + '</GBSet>')
                render(p/'tree', p, metadata, p/'figure')
                with (p/'figure.tips.tsv').open() as f:
                    tips = list(csv.DictReader(f, delimiter='\t'))
                self.assertEqual([tip['taxon_or_organism'] for tip in tips], list(taxa.values()))
                self.assertEqual(tips[0]['accession'], 'PROC-0' if source == 'bold' else 'REF0.1')
