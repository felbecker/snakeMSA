#!/usr/bin/env bash
# Remove failed alignments (see find_failed.sh) together with everything that
# depends on them, so that snakemake regenerates them on the next run: the
# score, benchmarks, intermediate files and the tool's results file, which
# concat_scores would otherwise consider up to date. Logs are kept, snakemake
# overwrites them on the rerun.
#
# Run from the snakeMSA root.
#
# Usage: tools/remove_failed.sh [-n] <benchmark name>
#   -n  dry run, only print what would be removed
set -euo pipefail

dry=0
if [[ "${1:-}" == "-n" ]]; then
    dry=1
    shift
fi
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 [-n] <benchmark name>" >&2
    exit 1
fi
name=$1

remove() {
    for path in "$@"; do
        if [[ -e "$path" || -L "$path" ]]; then
            echo "remove $path"
            if [[ $dry -eq 0 ]]; then
                rm -rf -- "$path"
            fi
        fi
    done
}

tools_dir=$(dirname "$0")
count=0
declare -A tools
while IFS= read -r msa; do
    [[ -n "$msa" ]] || continue
    # outputs/<name>/<tool>/alignments/<sample>
    sample=$(basename "$msa")
    tool_dir=$(dirname "$(dirname "$msa")")
    tools[$(basename "$tool_dir")]=1
    count=$((count + 1))
    remove "$msa" \
        "$tool_dir/scores/$sample.score" \
        "$tool_dir/benchmarks/$sample.txt" \
        "$tool_dir/benchmarks/$sample.gpu" \
        "$tool_dir/benchmarks/$sample.time" \
        "$tool_dir/projections/$sample" \
        "$tool_dir/workdirs/$sample" \
        "$tool_dir/trees/$sample.dnd"
done < <("$tools_dir/find_failed.sh" "$name")

for tool in "${!tools[@]}"; do
    remove "results/$name/$tool.out"
done

echo "$count failed alignment(s)$([[ $dry -eq 1 ]] && echo ' (dry run)')"
