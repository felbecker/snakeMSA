# Run inside the muscle env: conda run -n muscle bash run/muscle3D.sh
# reseek comes from bin/ (install.sh), the muscle env does not have it.
pdb="$1"
output="$2"
workdir="$3"
reseek="$(dirname "$0")/../bin/reseek"
"$reseek" -pdb2mega "$pdb" -output "$workdir/structs.mega" && \
muscle -align "$workdir/structs.mega" -output "$output" && \
python3 "$(dirname "$0")/../tools/fix_chain_ids.py" --msa "$output" --pdb "$pdb"
