# Run inside the muscle env: conda run -n muscle bash run/muscle3D_super7.sh
# reseek comes from bin/ (install.sh), as the env's 2.6.1 writes an empty
# distance matrix.
pdb="$1"
output="$2"
workdir="$3"
reseek="$(dirname "$0")/../bin/reseek"
# for up to ~10,000 structures
"$reseek" -convert "$pdb" -bca "$workdir/structs.bca" && \
"$reseek" -pdb2mega "$workdir/structs.bca" -output "$workdir/structs.mega" && \
"$reseek" -distmx "$workdir/structs.bca" -output "$workdir/structs.distmx" && \
muscle -super7 "$workdir/structs.mega" -distmxin "$workdir/structs.distmx" -reseek -output "$output" && \
python3 "$(dirname "$0")/../tools/fix_chain_ids.py" --msa "$output" --pdb "$pdb"
