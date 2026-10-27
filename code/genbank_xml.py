"""Read one or multiple concatenated GenBank XML documents from EDirect batches."""
import re
import xml.etree.ElementTree as ET


def iter_genbank_records(path):
    parser = ET.XMLPullParser(events=['end'])
    parser.feed('<EDirectBatches>')
    with open(path) as handle:
        for line in handle:
            # EDirect emits these declarations on single lines for every batch.
            line = re.sub(r'<\?xml\b.*?\?>', '', line)
            line = re.sub(r'<!DOCTYPE\s+GBSet\b[^>]*>', '', line)
            parser.feed(line)
            for _, element in parser.read_events():
                if element.tag == 'GBSeq':
                    yield element
                    element.clear()
                elif element.tag == 'GBSet':
                    element.clear()
    parser.feed('</EDirectBatches>')
    parser.close()
