#!/bin/bash
# Sets up the aligners:
# - a conda environment <name> for each envs/<name>.yaml
# - tools whose conda packages are outdated are installed to bin/ instead
#
# Usage: ./install.sh [name ...]
# Without names, everything is installed. Otherwise only the given conda
# environments and/or binaries (famsa, kalign, reseek).
set -euo pipefail

FAMSA_VERSION=2.5.2
KALIGN_VERSION=3.6.0
RESEEK_VERSION=2.8

repo_dir=$(cd "$(dirname "$0")" && pwd)
bin_dir="$repo_dir/bin"
targets=("$@")

known=(famsa kalign reseek)
for file in "$repo_dir"/envs/*.yaml; do
    known+=("$(basename "$file" .yaml)")
done
for target in "${targets[@]}"; do
    if [[ " ${known[*]} " != *" $target "* ]]; then
        echo "Unknown name $target. Choose from: ${known[*]}"
        exit 1
    fi
done

selected() {
    [ ${#targets[@]} -eq 0 ] || [[ " ${targets[*]} " == *" $1 "* ]]
}

create_env() {
    local file=$1 name
    name=$(basename "$file" .yaml)
    selected "$name" || return 0
    echo "Creating conda environment $name from $file"
    if conda info --envs | grep -q "^$name\s"; then
        echo "Environment $name already exists. Skipping."
    else
        conda env create --name "$name" --file "$file" || {
            echo "Failed to create environment $name. Exiting."
            exit 1
        }
    fi
}

# bioconda lags behind the GitHub releases, which ship static binaries.
install_famsa() {
    selected famsa || return 0
    if "$bin_dir/famsa" 2>&1 | grep -q "version $FAMSA_VERSION"; then
        echo "FAMSA $FAMSA_VERSION already installed. Skipping."
        return
    fi
    local platform
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) platform=x64_linux ;;
        Linux-aarch64) platform=arm64_linux ;;
        Darwin-x86_64) platform=x64_mac ;;
        Darwin-arm64) platform=arm64_mac ;;
        *) echo "No FAMSA binary for $(uname -s) $(uname -m). Exiting."; exit 1 ;;
    esac
    echo "Downloading FAMSA $FAMSA_VERSION ($platform) to $bin_dir"
    local tmp
    tmp=$(mktemp -d)
    curl -fsSL "https://github.com/refresh-bio/FAMSA/releases/download/v$FAMSA_VERSION/famsa-v$FAMSA_VERSION-$platform.tar.gz" \
        | tar -xz -C "$tmp"
    install -m 755 "$tmp/famsa" "$bin_dir/famsa"
    rm -rf "$tmp"
}

# bioconda lags behind GitHub, which has no binaries, so build from source.
install_kalign() {
    selected kalign || return 0
    if "$bin_dir/kalign" --version 2>&1 | grep -q "kalign $KALIGN_VERSION"; then
        echo "Kalign $KALIGN_VERSION already installed. Skipping."
        return
    fi
    local tool
    for tool in curl cmake make cc; do
        command -v "$tool" > /dev/null || {
            echo "Building Kalign requires $tool. Exiting."
            exit 1
        }
    done
    echo "Building Kalign $KALIGN_VERSION to $bin_dir"
    local tmp
    tmp=$(mktemp -d)
    curl -fsSL "https://github.com/TimoLassmann/kalign/archive/refs/tags/v$KALIGN_VERSION.tar.gz" \
        | tar -xz -C "$tmp"
    # static build, so that bin/kalign does not depend on a libkalign.so
    cmake -S "$tmp/kalign-$KALIGN_VERSION" -B "$tmp/build" \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF > "$tmp/build.log" 2>&1 \
        && cmake --build "$tmp/build" --target kalign-bin -j "$(nproc 2>/dev/null || echo 4)" \
            >> "$tmp/build.log" 2>&1 || {
        echo "Failed to build Kalign, see $tmp/build.log. Exiting."
        exit 1
    }
    install -m 755 "$tmp/build/src/kalign" "$bin_dir/kalign"
    rm -rf "$tmp"
}

# bioconda has reseek 2.6.1, whose -distmx writes no distances, so that
# muscle -super7 fails (rcedgar/muscle#85). 2.8 is the newest version that
# still has -pdb2mega and -distmx. The GitHub binaries need AVX-512, which
# Zen 2 lacks, so build from source for any AVX2 machine (x86-64-v3) instead
# of -march=native, which would tie the binary to the build host.
install_reseek() {
    selected reseek || return 0
    if "$bin_dir/reseek" -version 2>&1 | grep -q "reseek v$RESEEK_VERSION\."; then
        echo "reseek $RESEEK_VERSION already installed. Skipping."
        return
    fi
    local tool
    for tool in curl python3 make g++; do
        command -v "$tool" > /dev/null || {
            echo "Building reseek requires $tool. Exiting."
            exit 1
        }
    done
    local march=x86-64-v3
    [ "$(uname -m)" = x86_64 ] || march=native
    echo "Building reseek $RESEEK_VERSION ($march) to $bin_dir"
    local tmp src
    tmp=$(mktemp -d)
    curl -fsSL "https://github.com/rcedgar/reseek/archive/refs/tags/v$RESEEK_VERSION.tar.gz" \
        | tar -xz -C "$tmp"
    src="$tmp/reseek-$RESEEK_VERSION/src"
    # the build script generator, pinned by the release
    curl -fsSL "https://raw.githubusercontent.com/rcedgar/vcxproj_make/$(cat "$src/vcxproj_make_commit.txt")/vcxproj_make.py" \
        > "$tmp/vcxproj_make.py"
    echo '"0"' > "$src/gitver.txt"
    (cd "$src" \
        && python3 "$tmp/vcxproj_make.py" --openmp --lrt --pthread --nomake \
        && sed -i "s/-march=native/-march=$march/g" Makefile \
        && make -j "$(nproc 2>/dev/null || echo 4)") > "$tmp/build.log" 2>&1 || {
        echo "Failed to build reseek, see $tmp/build.log. Exiting."
        exit 1
    }
    install -m 755 "$tmp/reseek-$RESEEK_VERSION/bin/reseek" "$bin_dir/reseek"
    rm -rf "$tmp"
}

for file in "$repo_dir"/envs/*.yaml; do
    create_env "$file"
done

mkdir -p "$bin_dir"
install_famsa
install_kalign
install_reseek
