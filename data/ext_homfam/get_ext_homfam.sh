#!/bin/bash

if [ -d "unaligned" ]; then
    echo "unaligned/ already exists, skipping download and extraction."
else
    wget https://zenodo.org/records/15883613/files/extHomFam-v37.0.tar.gz
    tar -xzf extHomFam-v37.0.tar.gz
    mv extHomFam-v37.0/families unaligned/
    mv extHomFam-v37.0/references aligned/
fi

for file in unaligned/*; do
    python uniquify.py "$file"
done
for file in aligned/*; do
    python uniquify.py "$file"
done

# Final postprocessing: remove write permission for file owner
chmod u-w unaligned/*
chmod u-w aligned/*