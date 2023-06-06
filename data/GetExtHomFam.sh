#!/bin/bash
wget https://zenodo.org/record/6524237/files/extHomFam-v2.zip
unzip extHomFam-v2.zip
X=("small" "medium" "large" "xlarge" "huge")
for folder in "${X[@]}"
do
    mkdir ext_homfam_$folder
    mkdir ext_homfam_$folder/refs
    mkdir ext_homfam_$folder/train
    mv extHomFam-v2/$folder/* ext_homfam_$folder/train/
    for file in ext_homfam_$folder/train/*
    do
        filename=$(basename "$file")
        mv "extHomFam-v2/ref/$filename" "ext_homfam_$folder/refs"
    done
done
rm -r extHomFam-v2
rm extHomFam-v2.zip
python3 uniquify.py ext_homfam_large/train/DNA_pol3_beta
python3 uniquify.py ext_homfam_xlarge/train/ldh
python3 uniquify.py ext_homfam_huge/train/gtp
python3 uniquify.py ext_homfam_huge/train/thiored
python3 uniquify.py ext_homfam_huge/train/Epimerase
python3 uniquify.py ext_homfam_huge/train/apbact