# Run inside the foldmason env: conda run -n foldmason bash run/foldmason.sh
foldmason easy-msa "$1" "$2" "$3" && \
mv "${2}_aa.fa" "$2"
