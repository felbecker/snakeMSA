find outputs/*/*/alignments -type f -exec awk 'FNR==1 && /^FAILED\r?$/ { print FILENAME }' {} +
