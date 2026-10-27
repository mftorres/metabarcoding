"""Stream taxon-matched records from local gzipped BOLD FASTA databases."""
import csv
import gzip
from pathlib import Path
import re


def normalise(value):
    return re.sub(r'\s+', ' ', value.replace('_', ' ')).strip().casefold()


def database_files(path):
    path = Path(path)
    files = [path] if path.is_file() else sorted(p for p in path.glob('*') if p.name.endswith(('.fasta.gz', '.fa.gz', '.fna.gz')))
    if not files:
        raise FileNotFoundError(f'No gzipped BOLD FASTA database found at {path}')
    return files


def extract_bold(files, taxon, output, metadata, identifiers):
    needle = normalise(taxon)
    if not needle or taxon.strip().isdigit():
        raise ValueError('--bold requires a taxon name, not a numeric taxid')
    count = 0
    def records(path):
        header, sequence = None, []
        with gzip.open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('>'):
                    if header is not None:
                        yield header, ''.join(sequence)
                    header, sequence = line[1:].rstrip('\r\n'), []
                elif header is not None and needle in normalise(header):
                    sequence.append(line.strip())
            if header is not None:
                yield header, ''.join(sequence)
    # Stream gzip records; retain original BOLD metadata under safe tool-compatible IDs.
    with Path(output).open('x') as fasta, Path(metadata).open('x', newline='') as audit, Path(identifiers).open('x') as ids:
        writer = csv.writer(audit, delimiter='\t')
        writer.writerow(['sequence_id','process_id','marker','taxon','original_header','source_file'])
        for path in files:
            for header, sequence in records(path):
                if needle not in normalise(header):
                    continue
                if not sequence:
                    raise ValueError(f'Empty matching BOLD sequence: {header}')
                fields = header.split('|')
                if len(fields) < 4:
                    raise ValueError(f'Expected BOLD process|marker|country|taxonomy header: {header}')
                lineage = [v.strip() for v in fields[3].split(',') if v.strip().casefold() not in ('', 'none', 'nan', 'null', 'na')]
                label = lineage[-1] if lineage else 'unassigned'
                count += 1
                identifier = f'BOLD_{count:09d}'
                fasta.write(f'>{identifier} {header}\n')
                for start in range(0,len(sequence),80): fasta.write(sequence[start:start+80]+'\n')
                writer.writerow([identifier,fields[0],fields[1],label,header,str(Path(path).resolve())])
                ids.write(identifier+'\n')
    if not count:
        raise ValueError(f'No BOLD FASTA headers contain {taxon!r}')
    return count
