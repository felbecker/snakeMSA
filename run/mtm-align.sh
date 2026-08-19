
input_list="$3/mtm_align_input_list.txt"

find "$1" -maxdepth 1 -type f -name "*.pdb" | sort > "$input_list"

conda run -n mtm-align mtm-align -i "$input_list" -outdir "$3" && \
sed -i 's/^>\(.*\)\.pdb$/>\1/' "$3/result.fasta" && \
mv "$3/result.fasta" "$2"
