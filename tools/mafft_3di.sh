#!/usr/bin/env bash
# Runs MAFFT on 3Di sequences with foldseek's 3Di substitution matrix and
# outputs the alignment in amino acids.
# Usage: mafft_3di.sh [mafft options] input_3di.fasta input_aa.fasta > output.fasta
#
# tools/mat3di_mafft.out (from tools/make_mafft_3di_matrix.py) names the
# 3Di states with other letters so that MAFFT's amino acid groups for 6-mer
# distances hold similar 3Di states. The input is translated to these letters.
# FFT anchoring is turned off, as it uses amino acid polarity and volume.
# The gap pattern of the alignment is then put onto the amino acid sequences
# (same IDs and lengths as the 3Di ones), which the scoring matches to the
# structures.
set -euo pipefail

matrix="$(dirname "$0")/mat3di_mafft.out"
states=$(sed -n 's/^# 3Di letters: *//p' "$matrix")
letters=$(sed -n 's/^# MAFFT letters: *//p' "$matrix")
input="${@: -2:1}"
aa_input="${!#}"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
sed "/^>/!y/$states/$letters/" "$input" > "$tmp/input.fasta"
mafft "${@:1:$#-2}" --nofft --aamatrix "$matrix" "$tmp/input.fasta" \
    > "$tmp/aligned.fasta"
python3 "$(dirname "$0")/propagate_msa.py" "$tmp/aligned.fasta" "$aa_input"
