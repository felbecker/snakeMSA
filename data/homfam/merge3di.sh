#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MERGE_FASTA="$SCRIPT_DIR/../../3DUtil/merge_fasta.py"

PREDICTED_3DI="$SCRIPT_DIR/predicted_3di"
PDB_3DI="$SCRIPT_DIR/pdb_3di"
HOMSTRAD_3DI="$SCRIPT_DIR/../homstrad/3Di"
OUTPUT_DIR="$SCRIPT_DIR/3di"

mkdir -p "$OUTPUT_DIR"

for file in "$PREDICTED_3DI"/*.fasta; do
    family=$(basename "$file" .fasta)
    tmp=$(mktemp)

    # TODO: pdb_3di is currently broken because of unresolved residues
    # # Step 1: start from predicted_3di, override with pdb_3di
    # pdb_file="$PDB_3DI/$family.fasta"
    # if [ -f "$pdb_file" ]; then
    #     python3 "$MERGE_FASTA" -i "$file" --override_fasta "$pdb_file" -o "$tmp"
    # else
    #     cp "$file" "$tmp"
    # fi

    cp "$file" "$tmp"

    # Step 2: override result with homstrad/3Di
    homstrad_file="$HOMSTRAD_3DI/$family.fasta"
    if [ -f "$homstrad_file" ]; then
        python3 "$MERGE_FASTA" -i "$tmp" --override_fasta "$homstrad_file" -o "$OUTPUT_DIR/$family.fasta"
    else
        cp "$tmp" "$OUTPUT_DIR/$family.fasta"
    fi

    rm "$tmp"
    echo "Merged: $family"
done

echo "Done. Output in $OUTPUT_DIR"
