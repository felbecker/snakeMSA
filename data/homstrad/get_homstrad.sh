#!/bin/bash
HOMSTRAD_URL="https://homstrad.mizuguchilab.org/homstrad/data/"
HOMSRTAD_FILE="homstrad_with_PDB_2025_Dec_1.tar.gz"
DIR_NAME="homstrad_with_PDB"
UTIL_DIR="../../3DUtil"

if [ -f "$HOMSRTAD_FILE" ]; then
    echo "File $HOMSRTAD_FILE already exists. Skipping download."
else
    wget "$HOMSTRAD_URL$HOMSRTAD_FILE"
fi

if [ -d "$DIR_NAME" ]; then
    echo "Directory $DIR_NAME already exists. Skipping download."
else
    mkdir -p "$DIR_NAME"
    tar -xvzf "$HOMSRTAD_FILE" -C "$DIR_NAME"
fi

if [ -d aligned ] && [ -d pdbs ] && [ -d unaligned ]; then
    echo "Output directories already exist. Skipping creation."
else
    python3 clean_homstrad.py $DIR_NAME --min-sequences 5
fi

if [-d 3Di ]; then
    echo "Directory 3Di already exists. Skipping creation."
else
    mkdir -p 3Di
    for folder in pdbs/*; do
        python "$UTIL_DIR"/make_3di.py --pdb-dir "$folder" --output 3Di/$(basename "$folder").fasta
    done
fi
