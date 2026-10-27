#!/usr/bin/env bash
set -euo pipefail
# Post-classification only: never launches filtering, translation or classification.
# Activate taxonomy, then provide the output directory and optional Python CLI flags.
exec python code/run_sintax_notebooks.py "$@"
