pdb="$1"
output="$2"
workdir="$3"
# for up to ~10,000 structures
conda run -n muscle_reseek28 reseek -convert "$pdb" -bca "$workdir/structs.bca" && \
conda run -n muscle_reseek28 reseek -pdb2mega "$workdir/structs.bca" -output "$workdir/structs.mega" && \
conda run -n muscle_reseek28 reseek -distmx "$workdir/structs.bca" -output "$workdir/structs.distmx" && \
conda run -n muscle muscle -super7 "$workdir/structs.mega" -distmxin "$workdir/structs.distmx" -reseek -output "$output" && \
python3 tools/fix_chain_ids.py --msa "$output" --pdb "$pdb"