#!/bin/bash
# Find PDB structures and extract 3Di sequences for all homfam families.
# Run from the snakeMSA root directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UNALIGNED="$SCRIPT_DIR/unaligned"
ALIGNED="$SCRIPT_DIR/aligned"
OUTPUT_DIR="$SCRIPT_DIR/pdb_3di"
WORKDIR="$SCRIPT_DIR/workdir"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORKDIR"

UNALIGNED_STATS="$SCRIPT_DIR/statistics_unaligned.tsv"
ALIGNED_STATS="$SCRIPT_DIR/statistics_aligned.tsv"

for file in "$UNALIGNED"/*.vie; do
    stem="$(basename "${file%.vie}")"
    output="$OUTPUT_DIR/${stem}.fasta"
    if [ -f "$output" ]; then
        echo "Skipping $file (output already exists)"
        continue
    fi
    echo "Processing $file"
    python 3DUtil/find_structures.py \
        --input "$file" \
        --workdir "$WORKDIR" \
        --3di-output "$output"

    python 3DUtil/pdb_statistics.py --input "$output" --reference "$file" --append-tsv "$UNALIGNED_STATS"
    if [ -f "$ALIGNED/${stem}.vie" ]; then
        python 3DUtil/pdb_statistics.py --input "$output" --reference "$ALIGNED/${stem}.vie" --append-tsv "$ALIGNED_STATS"
    fi
done
