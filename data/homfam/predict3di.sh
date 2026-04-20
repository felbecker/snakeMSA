#!/bin/bash
# compute 3di predictions for the homfam dataset
util_dir=../../3DUtil
mkdir -p predicted_3di predicted_3di_esm
for file in unaligned/*.vie; do
    base="$(basename "${file%.vie}.fasta")"

    output_prostt5="predicted_3di/$base"
    if [ -f "$output_prostt5" ]; then
        echo "Skipping prostt5 $file (output already exists)"
    else
        echo "Processing prostt5 $file"
        python "$util_dir/make_3di.py" \
        --fasta "$file" \
        --output "$output_prostt5" \
        --pred-model prostt5
    fi

    output_esm="predicted_3di_esm/$base"
    if [ -f "$output_esm" ]; then
        echo "Skipping esm2 $file (output already exists)"
    else
        echo "Processing esm2 $file"
        python "$util_dir/make_3di.py" \
        --fasta "$file" \
        --output "$output_esm" \
        --pred-model esm2
    fi
done