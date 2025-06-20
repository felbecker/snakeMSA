find outputs/$1/*/alignments -type f -exec awk 'FNR==1 && /^FAILED\r?$/ { print FILENAME }' {} +
