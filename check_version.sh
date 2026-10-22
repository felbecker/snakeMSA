#!/bin/bash
# Prints the versions of the tools installed by install.sh and where they are.
#
# Usage: ./check_version.sh
set -uo pipefail

repo_dir=$(cd "$(dirname "$0")" && pwd)
bin_dir="$repo_dir/bin"
envs=$(conda info --envs)

report() {
    printf "%-12s %-28s %s\n" "$1" "${2:-not installed}" "$3"
}

# Runs a version command inside a conda environment and keeps the first
# non-empty line of its output.
env_version() {
    local env=$1 cmd=$2
    grep -q "^$env\s" <<< "$envs" || return 0
    conda run -n "$env" bash -c "$cmd" 2>&1 | grep -v '^\s*$' | head -n 1
}

bin_version() {
    local tool=$1 cmd=$2
    [ -x "$bin_dir/$tool" ] || return 0
    bash -c "$cmd" 2>&1 | grep -v '^\s*$' | head -n 1
}

report tool version location
report famsa "$(bin_version famsa "'$bin_dir/famsa' 2>&1 | sed -n 's/^ *version //p'")" \
    "bin/famsa (GitHub release)"
report kalign "$(bin_version kalign "'$bin_dir/kalign' --version | sed 's/^kalign //'")" \
    "bin/kalign (built from source)"
report reseek "$(bin_version reseek "'$bin_dir/reseek' -version | sed 's/^reseek v//; s/\\.i86.*//'")" \
    "bin/reseek (built from source)"
report clustalo "$(env_version clustalo "clustalo --version")" "conda env clustalo"
report mafft "$(env_version mafft "mafft --version")" "conda env mafft"
report muscle "$(env_version muscle "muscle -version | cut -d ' ' -f 2")" "conda env muscle"
report t-coffee "$(env_version t-coffee \
    "t_coffee -version | sed -n 's/.*Version_\([^ ]*\).*/\1/p'")" "conda env t-coffee"
report learnMSA "$(env_version learnMSA \
    "\$CONDA_PREFIX/bin/python -c 'import importlib.metadata as m; print(m.version(\"learnMSA\"))'")" \
    "conda env learnMSA"
report snakemake "$(env_version snakeMSA "snakemake --version")" "conda env snakeMSA"
report foldmason "$(env_version foldmason "foldmason version")" "conda env foldmason"
