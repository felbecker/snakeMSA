#!/bin/bash
# learnMSA on per-residue ProstT5 3Di probabilities.
#
# Usage: run/run_prostT5_learnMSA.sh INPUT OUTPUT WORKDIR [learnMSA args...]
#
# 1. learnMSA-3di predicts per-residue 3Di logits for INPUT with ProstT5 and
#    writes WORKDIR/<name>.3di.npz (plus the argmax letters as
#    WORKDIR/<name>.3di.fasta). An existing npz is reused, so other runs can
#    point --struct at it.
# 2. learnMSA aligns INPUT with --struct WORKDIR/<name>.3di.npz and any extra
#    arguments.
#
# Both steps use the learnMSA worktree ~/src/learnMSA-3di-uncertainty
# (branch 3di_uncertainty) in the conda env learnMSAtorch.
set -euo pipefail

input="$1"
output="$2"
workdir="$3"
shift 3

learnmsa="$(cd "$(dirname "$0")/../../learnMSA-3di-uncertainty" && pwd)"
name="$(basename "${input%.*}")"
npz="$workdir/$name.3di.npz"
mkdir -p "$workdir"

if [ -f "$npz" ]; then
    echo "Reusing 3Di logits $npz"
else
    # Write under a temporary name first so that an interrupted run never
    # leaves a truncated npz behind for later runs to reuse.
    tmp="$workdir/$name.3di.tmp.npz"
    PYTHONPATH="$learnmsa" conda run --live-stream -n learnMSAtorch \
        python -m learnMSA.structure.predict_3di \
        -i "$input" -o "$tmp" --fasta "$workdir/$name.3di.fasta"
    mv "$tmp" "$npz"
fi

conda run --live-stream -n learnMSAtorch python "$learnmsa/learnMSA.py" \
    -i "$input" -o "$output" --struct "$npz" --work_dir "$workdir" "$@"
