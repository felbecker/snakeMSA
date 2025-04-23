#!/bin/bash
wget http://www.clustal.org/omega/homfam-20110613-25.tar.gz
tar -xvzf homfam-20110613-25.tar.gz
mkdir aligned
mkdir unaligned
# move all files that match *_ref.vie into aligned and strip the _ref 
# move the rest with the pattern *_test-only.vie into unaligned and strip "_test-only"
for file in *_ref.vie; do
    mv "$file" "aligned/$(basename "$file" _ref.vie).vie"
done
for file in *_test-only.vie; do
    mv "$file" "unaligned/$(basename "$file" _test-only.vie).vie"
done
rm -r homfam-20110613-25.tar.gz