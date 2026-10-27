import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET
from genbank_xml import iter_genbank_records
from plot_phylogeny import reference_names
from postprocess_sintax import validate_outputs

class BatchedXMLTests(unittest.TestCase):
    def test_multiple_documents(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'reference_records.gb.xml'
            path.write_text(''.join('<?xml version="1.0"?>\n<!DOCTYPE GBSet PUBLIC "x" "x">\n<GBSet><GBSeq>'
                f'<GBSeq_accession-version>{acc}</GBSeq_accession-version><GBSeq_organism>Tipula {acc}</GBSeq_organism>'
                '</GBSeq></GBSet>\n' for acc in ('ONE.1','TWO.1')))
            validate_outputs([path])
            self.assertEqual(len(list(iter_genbank_records(path))),2)
            self.assertEqual(reference_names(path),{'ONE.1':'Tipula ONE.1','TWO.1':'Tipula TWO.1'})
            path.write_text(path.read_text()+'<GBSet><GBSeq>')
            with self.assertRaises(ET.ParseError): list(iter_genbank_records(path))
