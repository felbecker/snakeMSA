#!/bin/bash
# compute 3di predictions for the homfam dataset
mkdir -p predicted_3di
for file in unaligned/*.vie; do
    output="predicted_3di/$(basename "${file%.vie}.fasta")"
    if [ -f "$output" ]; then
        echo "Skipping $file (output already exists)"
        continue
    fi
    echo "Processing $file"
    python 3DUtil/make_3di.py \
    --fasta "$file" \
    --output "$output" \
    --prostt5-model "3DUtil/prostT5"
done