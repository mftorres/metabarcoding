#!/usr/bin/env bash
set -euo pipefail
# Convenience entry point; pass any sample folders as arguments.
exec bash code/run_metabarcoding_pipeline.sh "$@"
