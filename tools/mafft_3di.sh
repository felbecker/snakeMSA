#!/usr/bin/env bash
# Runs MAFFT on 3Di sequences with foldseek's 3Di substitution matrix.
# Usage: mafft_3di.sh [mafft options] input.fasta > output.fasta
#
# tools/mat3di_mafft.out (from tools/make_mafft_3di_matrix.py) names the
# 3Di states with other letters so that MAFFT's amino acid groups for 6-mer
# distances hold similar 3Di states. The input is translated to these letters
# and the alignment back. FFT anchoring is turned off, as it uses amino acid
# polarity and volume.
set -euo pipefail

matrix="$(dirname "$0")/mat3di_mafft.out"
states=$(sed -n 's/^# 3Di letters: *//p' "$matrix")
letters=$(sed -n 's/^# MAFFT letters: *//p' "$matrix")
input="${!#}"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
sed "/^>/!y/$states/$letters/" "$input" > "$tmp/input.fasta"
mafft "${@:1:$#-1}" --nofft --aamatrix "$matrix" "$tmp/input.fasta" \
    | sed "/^>/!y/$letters/$states/"
