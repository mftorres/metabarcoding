"""Local CLI adapter for the installed kc-align mixed-mode thread forwarding bug."""
import argparse
import kcalign
from kcalign_parallel import create_lists


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['mixed'], required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--sequences', required=True)
    parser.add_argument('--table', type=int, default=5)
    parser.add_argument('--keep', action='store_true')
    parser.add_argument('--threads', type=int, default=6)
    args = parser.parse_args()
    if args.threads < 1 or (args.threads != 1 and args.threads % 3):
        parser.error('Mixed mode requires 1 thread or a positive multiple of 3')
    # Correct the installed parallel rejection bug without modifying the conda package.
    kcalign.create_lists = create_lists
    return kcalign.mixed_mode(args.reference, args.sequences, False, args.table, args.keep, 'none', args.threads)

if __name__ == '__main__':
    raise SystemExit(main())
