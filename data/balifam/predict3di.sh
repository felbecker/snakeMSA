#!/bin/bash
# compute 3di predictions for the balifam dataset
script_dir="$(cd "$(dirname "$0")" && pwd)"
util_dir="$script_dir/../../3DUtil"
for dataset_dir in balifam/balifam100 balifam/balifam1000 balifam/balifam10000; do
    pushd "$script_dir/$dataset_dir" > /dev/null
    mkdir -p predicted_3di predicted_3di_esm
    for file in unaligned/*; do
        base="$(basename "${file}.fasta")"

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
    done
    popd > /dev/null
done