#!/bin/bash
wget https://zenodo.org/record/6524237/files/extHomFam-v2.zip
python3 -m zipfile -e extHomFam-v2.zip .
X=("small" "medium" "large" "xlarge" "huge")
for folder in "${X[@]}"
do
    mkdir ext_homfam_$folder
    mkdir ext_homfam_$folder/aligned
    mkdir ext_homfam_$folder/unaligned
    mv extHomFam-v2/$folder/* ext_homfam_$folder/unaligned/
    for file in ext_homfam_$folder/unaligned/*
    do
        filename=$(basename "$file")
        mv "extHomFam-v2/ref/$filename" "ext_homfam_$folder/aligned"
    done
done
rm -r extHomFam-v2
rm extHomFam-v2.zip
python3 uniquify.py ext_homfam_large/unaligned/DNA_pol3_beta
python3 uniquify.py ext_homfam_xlarge/unaligned/ldh
python3 uniquify.py ext_homfam_huge/unaligned/gtp
python3 uniquify.py ext_homfam_huge/unaligned/thiored
python3 uniquify.py ext_homfam_huge/unaligned/Epimerase
python3 uniquify.py ext_homfam_huge/unaligned/apbact

# Final postprocessing: remove write permission for file owner
for folder in "${X[@]}"
do
    chmod u-w ext_homfam_$folder/unaligned/*
    chmod u-w ext_homfam_$folder/aligned/*
done
