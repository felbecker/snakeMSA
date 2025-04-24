#!/bin/bash
wget http://www.clustal.org/omega/homfam-20110613-25.tar.gz
tar -xvzf homfam-20110613-25.tar.gz
mkdir aligned
mkdir unaligned
# move all files that match *_ref.vie into aligned and strip the _ref 
# move the rest with the pattern *_test-only.vie into unaligned and strip "_test-only"
# also append the reference sequences to the test sequences and remove gaps
for file in *_ref.vie; do
    mv "$file" "aligned/$(basename "$file" _ref.vie).vie"
done
for file in *_test-only.vie; do
    # add references to the test sequences
    target_file=$(basename "$file" _test-only.vie).vie
    mv "$file" "unaligned/$target_file"
    cat "unaligned/$target_file" "aligned/$target_file" > "unaligned/$target_file.tmp"
    # remove gaps and keep file permissions
    sed -i '/^>/!s/[-.]//g' "unaligned/$target_file.tmp" #remove gaps
    # Save original permissions
    PERMS=$(stat -c %a unaligned/$target_file)
    # Make the file writable (temporarily)
    chmod u+w unaligned/$target_file
    # Replace it
    mv -f "unaligned/$target_file.tmp" "unaligned/$target_file"
    # Restore original permissions
    chmod $PERMS "unaligned/$target_file"
done
rm -r homfam-20110613-25.tar.gz