"""Isolated migration checks: temporary fixtures and fake external commands only."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from sintax_parser import discover_outputs, parse_output, load_outputs
from midori_taxonomy import catalog

CODE = Path(__file__).resolve().parent


class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='sintax_test_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.sample = self.root / 'input' / 'unclassified'
        self.sample.mkdir(parents=True)
        self.output_root = self.root / 'results'
        self.output = self.output_root / 'unclassified' / 'unclassified_sintaxSP.tsv'
        self.output.parent.mkdir(parents=True)

    def test_parse_confidence_provenance_and_empty_calls(self):
        self.output.write_text(
            'split-by-adapter-left-r1 parent_read_id=r1\tk:Eukaryota_2759(1.00),c:Insecta_50557(0.80),s:Fly one_123(0.30)\t+\tk:Eukaryota_2759\n'
            'r2\t\t+\t\n'
            'r3\tc:Insecta_50557(1.00)\t-\n')
        rows = list(parse_output(self.output))
        self.assertEqual(rows[0]['raw_ont_read_id'], 'r1')
        self.assertEqual(rows[0]['sintax_taxid'], '2759')
        self.assertEqual(rows[0]['accepted_class'], '')
        self.assertEqual(rows[0]['class'], 'Insecta_50557')
        self.assertEqual(rows[0]['source_folder'], 'unclassified')
        self.assertEqual(rows[0]['source_name'], self.output.name)
        self.assertFalse(rows[1]['classified'])
        self.assertFalse(rows[2]['classified'])
        frame = load_outputs([self.output])
        self.assertEqual(str(frame['sintax_taxid'].dtype), 'Int64')
        (self.output.parent / 'unclassified_consensus_sintaxSP.tsv').write_text('ignored')
        (self.output.parent / 'Tipulidae_sintaxSP.tsv').write_text('ignored')
        self.assertEqual(discover_outputs([self.output_root, self.output]), [self.output])
        other = self.root / 'other' / 'unclassified' / self.output.name
        other.parent.mkdir(parents=True)
        other.write_text(self.output.read_text())
        self.assertEqual(len(load_outputs([self.output, other])), 6)

    def test_malformed_input(self):
        for text in ['r\tbad\t+\n', 'r\tk:Eukaryota(2)\t+\n', 'r\n']:
            self.output.write_text(text)
            with self.assertRaises(ValueError):
                list(parse_output(self.output))

    def test_tracking_split_join_and_no_overwrite(self):
        raw = self.sample / 'reads.fastq.gz'
        with gzip.open(raw, 'wt') as handle:
            handle.write('@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n')
        self.output.write_text(
            'split-by-adapter-left-r1\tc:Insecta_50557(1.00)\t+\tc:Insecta_50557\n'
            'split-by-adapter-right-r1\tc:Insecta_50557(0.20)\t-\t\n')
        primers = self.root / 'primers.tsv'
        primers.write_text('merged_folder\tbarcode_folder\tseqID\tforward_patternName\ninput\tunclassified\tr1\tprimer1\nother\tunclassified\tr2\twrong\n')
        mft = self.root / 'mft.fasta'
        pcr = self.root / 'pcr.fasta'
        mft.write_text('>primer1\nACGT\n')
        pcr.write_text('>bc1001\nACGT\n')
        out = self.root / 'tracking' / 'unclassified_read_sintax.tsv'
        cmd = [sys.executable, str(CODE / 'build_barcode_read_tracking_table.py'),
               '--barcode-dir', str(self.sample), '--primer-best-hits', str(primers),
               '--sintax-output', str(self.output), '--output', str(out),
               '--mft-primers', str(mft), '--pcr-primers', str(pcr)]
        subprocess.run(cmd, check=True, capture_output=True)
        with out.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]['forward_PCR_primer'], 'bc1001')
        self.assertEqual(rows[-1]['forward_PCR_primer'], '')
        self.assertEqual(rows[-1]['sintax_taxid'], '')
        insecta = out.parent / 'unclassified_insecta_read_sintax_reduced.tsv'
        with insecta.open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter='\t'))), 1)
        original = out.read_bytes()
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
        self.assertEqual(out.read_bytes(), original)

    def setup_pipeline(self):
        binary = self.root / 'bin'
        binary.mkdir()
        fake = binary / 'vsearch'
        fake.write_text('''#!/usr/bin/env python3
import os,sys
from pathlib import Path
args=sys.argv[1:]
with open(os.environ['CALL_LOG'],'a') as f: f.write(' '.join(args)+'\\n')
if os.environ.get('FAIL_CLASSIFIER') == '1' and '--sintax' in args: sys.exit(9)
for flag in ('--output','--tabbedout'):
 if flag in args: Path(args[args.index(flag)+1]).write_text('fixture-output\\n')
''')
        fake.chmod(0o755)
        filtered = self.sample / 'filtered'
        filtered.mkdir()
        (filtered / 'asinput.fastq.gz').write_bytes(b'fixture, never read by real software')
        env = dict(os.environ, PATH=str(binary)+os.pathsep+os.environ['PATH'],
                   CONFIG_FILE=str(self.root / 'missing.env'), SINTAX_OUTPUT_ROOT=str(self.output_root),
                   SINTAX_REFERENCE=str(self.root / 'ref.fasta'), SINTAX_DB=str(self.root / 'ref.udb'),
                   SINTAX_SUFFIX='sintaxSP', CALL_LOG=str(self.root / 'calls'), FORCE='0', FORCE_SINTAX='0',
                   MAX_JOBS='2')
        return env

    def run_pipeline(self, env, *args):
        return subprocess.run(['bash', str(CODE / 'run_metabarcoding_pipeline.sh'), *args],
                              cwd=self.root, env=env, capture_output=True, text=True)

    def test_pipeline_skip_dryrun_and_database_reuse(self):
        env = self.setup_pipeline()
        self.output.write_text('existing result')
        before = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(self.run_pipeline(env, str(self.sample)).returncode, 0)
        self.assertEqual(self.run_pipeline(env, '--dry-run', str(self.sample)).returncode, 0)
        after = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before, after)
        self.output.unlink()  # temporary fixture only
        ref = Path(env['SINTAX_REFERENCE'])
        with gzip.open(str(ref)+'.gz', 'wt') as handle: handle.write('>ref\nACGT\n')
        result = self.run_pipeline(env, str(self.sample))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(ref.read_text(), '>ref\nACGT\n')
        calls = Path(env['CALL_LOG']).read_text()
        self.assertIn('--makeudb_usearch', calls)
        for flag in ['--sintax_cutoff 1', '--threads 8', '--randseed 666', '--strand both',
                     '--sintax_random', '--maxseqlength 2000', '--minseqlength 300']:
            self.assertIn(flag, calls)
        ref_time, db_time = ref.stat().st_mtime_ns, Path(env['SINTAX_DB']).stat().st_mtime_ns
        env['FORCE_SINTAX'] = '1'
        self.assertEqual(self.run_pipeline(env, str(self.sample)).returncode, 0)
        self.assertEqual(ref.stat().st_mtime_ns, ref_time)
        self.assertEqual(Path(env['SINTAX_DB']).stat().st_mtime_ns, db_time)
        self.assertEqual(Path(env['CALL_LOG']).read_text().count('--makeudb_usearch'), 1)
        env['FAIL_CLASSIFIER'] = '1'
        before = self.output.read_bytes()
        self.assertNotEqual(self.run_pipeline(env, str(self.sample)).returncode, 0)
        self.assertEqual(self.output.read_bytes(), before)
        self.assertFalse(list(self.output.parent.glob('*.tsv.*')))

    def test_existing_fasta_not_decompressed_and_duplicate_names(self):
        env = self.setup_pipeline()
        ref = Path(env['SINTAX_REFERENCE'])
        ref.write_text('>existing\nACGT\n')
        Path(str(ref)+'.gz').write_bytes(b'not gzip; must never be read')
        self.assertEqual(self.run_pipeline(env, str(self.sample)).returncode, 0)
        self.assertEqual(ref.read_text(), '>existing\nACGT\n')
        other = self.root / 'other' / 'unclassified'
        other.mkdir(parents=True)
        self.assertNotEqual(self.run_pipeline(env, str(self.sample), str(other)).returncode, 0)

    def test_midori_repeated_taxid_has_no_self_cycle(self):
        ref = self.root / 'ref.fasta'
        ref.write_text('>ref;tax=k:Eukaryota_2759,f:Flies_10,g:genus_Fly sp._12,s:Fly sp._12;\nACGT\n')
        names, parents, ranks = catalog(str(ref))
        self.assertEqual(names[12], 'Fly sp.')
        self.assertEqual(ranks[12], 'species')
        self.assertEqual(parents[12], 10)


if __name__ == '__main__':
    unittest.main()
