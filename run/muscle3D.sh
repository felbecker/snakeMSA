pdb="$1"
output="$2"
workdir="$3"
conda run -n muscle reseek -pdb2mega "$pdb" -output "$workdir/structs.mega" && \
conda run -n muscle muscle -align "$workdir/structs.mega" -output "$output" && \
python3 tools/fix_chain_ids.py --msa "$output" --pdb "$pdb"