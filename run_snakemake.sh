#!/usr/bin/env zsh
set -euo pipefail

#export CONDA_SUBDIR=osx-64

snakemake --use-conda "$@"
