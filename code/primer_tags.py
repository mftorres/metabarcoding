#!/usr/bin/env python3
"""Identify the forward and reverse primer tag of every read from a cutadapt info file.

cutadapt finds and trims the untagged primer cores. This step scores the tag block
next to each core (the 5-nt GGTAG spacer plus the 16-nt tag) against every known tag,
allowing substitutions, insertions and deletions. A tag is assigned only when its
best match has at most --max-tag-errors errors and beats the second-best tag by at
least --min-margin errors; otherwise the read is reported as ambiguous or too noisy.
"""
import argparse
import csv
import gzip
import sys

import numpy as np

from seqio import fasta

STATUSES = ('assigned', 'ambiguous', 'too_many_errors', 'truncated', 'primer_not_found')
COLUMNS = ['read_id', 'reverse_complemented'] + [
    f'{direction}_{field}' for direction in ('forward', 'reverse')
    for field in ('primer', 'status', 'core_errors', 'tag_errors', 'tag_margin', 'tag_candidates')]
CODES = np.full(256, 4, dtype=np.uint8)
for _code, _base in enumerate('ACGT'):
    CODES[ord(_base)] = CODES[ord(_base.lower())] = _code
PAD = 5  # flank bases absent beyond the searched window: aligning a tag base costs one error
OFF_READ = 6  # the read itself starts here: tag bases beyond the read end cost nothing
COMPLEMENT = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')


def load_primers(path, block_length):
    """Split each tagged primer into its tag block and shared core.

    Names ending in 'rc' are reverse primers. Every primer of one direction must
    share an identical core, and tag blocks must contain only A/C/G/T.
    """
    groups = {'forward': {}, 'reverse': {}}
    for header, sequence in fasta(path):
        name, sequence = header.split()[0], sequence.upper()
        direction = 'reverse' if name.endswith('rc') else 'forward'
        if name in groups[direction]:
            raise ValueError(f'Duplicate primer name in {path}: {name}')
        groups[direction][name] = sequence
    result = {}
    for direction, primers in groups.items():
        if len(primers) < 2:
            raise ValueError(f'{path}: need at least two {direction} tagged primers')
        cores = {sequence[block_length:] for sequence in primers.values()}
        if len(cores) != 1:
            raise ValueError(f'{path}: {direction} primers do not share one core after {block_length} nt')
        blocks = {name: sequence[:block_length] for name, sequence in primers.items()}
        if any(set(block) - set('ACGT') for block in blocks.values()):
            raise ValueError(f'{path}: {direction} tag blocks must contain only A/C/G/T')
        result[direction] = (cores.pop(), blocks)
    return result


def cutadapt_adapter(cores, core_errors):
    """Linked adapter: forward core at the 5' end, reverse-complemented reverse core at the 3' end."""
    forward, reverse = cores['forward'][0], cores['reverse'][0].translate(COMPLEMENT)[::-1]
    return (f'{forward};e={core_errors};min_overlap={len(forward)};optional...'
            f'{reverse};e={core_errors};min_overlap={len(reverse)};optional')


def encode(windows, width, tag_length=None):
    """Right-align windows in a uint8 matrix.

    Leading positions are PAD, or OFF_READ when the read starts inside the tag
    (window shorter than tag_length), so only the observed part of a tag is scored.
    """
    matrix = np.full((len(windows), width), PAD, dtype=np.uint8)
    for row, window in enumerate(windows):
        window = window[-width:]
        if tag_length is not None and len(window) < tag_length:
            matrix[row] = OFF_READ
        if window:
            matrix[row, width - len(window):] = CODES[np.frombuffer(window.encode(), dtype=np.uint8)]
    return matrix


def tag_distances(windows, tags):
    """Edit distance of each whole tag ending at the window end, with a free window start."""
    count, width = windows.shape
    length = tags.shape[1]
    previous = np.broadcast_to(np.arange(length + 1, dtype=np.int16)[:, None, None],
                               (length + 1, count, tags.shape[0])).copy()
    current = np.empty_like(previous)
    for column in range(width):
        current[0] = 0
        base = windows[:, column][:, None]
        for row in range(1, length + 1):
            cost = (tags[:, row - 1][None, :] != base) & (base != OFF_READ)
            np.add(previous[row - 1], cost, out=current[row])
            np.minimum(current[row], previous[row] + 1, out=current[row])
            np.minimum(current[row], current[row - 1] + 1, out=current[row])
        previous, current = current, previous
    return previous[length]


class TagScorer:
    """Score windows in batches; identical windows are scored once."""

    def __init__(self, blocks, reverse, max_errors, margin, slack, min_observed):
        self.names = list(blocks)
        # Reverse tags follow the core in the oriented read: align reversed sequences.
        sequences = [block.translate(COMPLEMENT) if reverse else block for block in blocks.values()]
        self.tags = encode(sequences, len(sequences[0]))
        self.width = self.tags.shape[1] + slack
        self.reverse, self.max_errors, self.margin = reverse, max_errors, margin
        self.min_observed = min_observed
        self.cache = {}

    def window(self, flank):
        return flank[:self.width][::-1] if self.reverse else flank[-self.width:]

    def score(self, flanks):
        windows = [self.window(flank) for flank in flanks]
        new = [w for w in dict.fromkeys(windows) if w not in self.cache]
        if len(self.cache) > 2_000_000:
            self.cache.clear()
        for start in range(0, len(new), 20_000):
            batch = new[start:start + 20_000]
            distances = tag_distances(encode(batch, self.width, self.tags.shape[1]), self.tags)
            order = np.argsort(distances, axis=1, kind='stable')
            for window, row, ranks in zip(batch, distances, order):
                best, second = int(row[ranks[0]]), int(row[ranks[1]])
                if len(window) < self.min_observed:
                    status, primer, candidates = 'truncated', '', ''
                elif best > self.max_errors:
                    status, primer, candidates = 'too_many_errors', '', ''
                elif second - best < self.margin:
                    status, primer = 'ambiguous', ''
                    candidates = ';'.join(self.names[i] for i in ranks if row[i] < best + self.margin)
                else:
                    status, primer, candidates = 'assigned', self.names[ranks[0]], ''
                self.cache[window] = (primer, status, best, second - best, candidates)
        return [self.cache[w] for w in windows]


def info_reads(handle):
    """Group consecutive cutadapt info lines per read: {'5': (errors, flank), '3': ...}."""
    current, parts = None, {}
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        read_id = fields[0].split()[0]
        if read_id != current:
            if current is not None:
                yield current, parts
            current, parts = read_id, {}
        if len(fields) < 12 or fields[1] == '-1':
            continue
        side = '5' if fields[7].endswith(';1') else '3'
        if side in parts:
            raise ValueError(f'Repeated {side} primer match for {read_id}')
        parts[side] = (int(fields[1]), fields[4] if side == '5' else fields[6], fields[11])
    if current is not None:
        yield current, parts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--primers', required=True, help='Tagged primer FASTA (bcXXXX forward, bcXXXXrc reverse)')
    parser.add_argument('--info', help='cutadapt --info-file to read (a file or named pipe)')
    parser.add_argument('--output', help='Per-read tag table (.tsv.gz)')
    parser.add_argument('--block-length', type=int, default=21, help='Spacer + tag length at the start of each primer')
    parser.add_argument('--max-tag-errors', type=int, default=4)
    parser.add_argument('--min-margin', type=int, default=2)
    parser.add_argument('--slack', type=int, default=4, help='Extra flank bases searched for tags with insertions')
    parser.add_argument('--min-observed', type=int, default=21,
                        help='Minimum tag-block bases present before the read starts; shorter flanks are "truncated"')
    parser.add_argument('--print-adapter', type=int, metavar='CORE_ERRORS',
                        help='Print the cutadapt linked-adapter specification and exit')
    args = parser.parse_args()
    primers = load_primers(args.primers, args.block_length)
    if args.print_adapter is not None:
        print(cutadapt_adapter(primers, args.print_adapter))
        return
    if not args.info or not args.output:
        parser.error('--info and --output are required')
    scorers = {side: TagScorer(primers[direction][1], side == '3', args.max_tag_errors, args.min_margin,
                               args.slack, args.min_observed)
               for side, direction in (('5', 'forward'), ('3', 'reverse'))}
    counts = {(direction, status): 0 for direction in ('forward', 'reverse') for status in STATUSES}
    seen = set()

    def flush(pending, writer):
        scores = {side: iter(scorer.score([parts[side][1] for _, parts in pending if side in parts]))
                  for side, scorer in scorers.items()}
        for read_id, parts in pending:
            flags = {value[2] for value in parts.values()}
            row = [read_id, 'yes' if '1' in flags else 'no' if flags else '']
            for side, direction in (('5', 'forward'), ('3', 'reverse')):
                if side in parts:
                    primer, status, errors, margin, candidates = next(scores[side])
                    row += [primer, status, parts[side][0], errors, margin, candidates]
                else:
                    status = 'primer_not_found'
                    row += ['', status, '', '', '', '']
                counts[direction, status] += 1
            writer.writerow(row)

    with open(args.info) as handle, gzip.open(args.output, 'wt', newline='') as output:
        writer = csv.writer(output, delimiter='\t', lineterminator='\n')
        writer.writerow(COLUMNS)
        pending = []
        for read_id, parts in info_reads(handle):
            if read_id in seen:
                raise ValueError(f'Read {read_id} appears twice in the cutadapt info file')
            seen.add(read_id)
            pending.append((read_id, parts))
            if len(pending) >= 50_000:
                flush(pending, writer)
                pending = []
        flush(pending, writer)
    total = len(seen)
    for direction in ('forward', 'reverse'):
        summary = ', '.join(f'{status} {counts[direction, status]:,} ({counts[direction, status] / max(total, 1):.1%})'
                            for status in STATUSES)
        print(f'{direction} tags in {total:,} reads: {summary}', file=sys.stderr)


if __name__ == '__main__':
    main()
