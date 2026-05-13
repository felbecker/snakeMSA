#!/bin/bash
BALIFAM_URL="https://github.com/rcedgar/balifam"
if [ -d "balifam" ]; then
    echo "Directory balifam already exists. Skipping clone."
else
    git clone "$BALIFAM_URL"
fi
cd balifam
for dir in balifam*; do
    mv "$dir/in" "$dir/unaligned"
    mv "$dir/ref" "$dir/aligned"
done
# Protect benchmark files from writing
chmod -w balifam*/unaligned/* balifam*/aligned/*
