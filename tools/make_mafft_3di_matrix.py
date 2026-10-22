"""Convert foldseek's 3Di substitution matrix into a MAFFT --aamatrix file.

MAFFT treats the 20 letters of a user matrix as amino acids beyond the
scores themselves:

1. Without a "frequency" line it normalizes the matrix (zero expected
   score, average self score of 600 that the gap penalties refer to) with
   BLOSUM62 amino acid frequencies. The output gets foldseek's 3Di
   background instead.
2. The 6-mer distances of the first guide tree reduce letters to the 6
   Dayhoff groups AGPST, ILMV, NDQE, RHK, FWY and C. The 3Di states are
   relabeled so that each group holds 3Di states that are often aligned to
   each other, i.e. the partition maximizes the probability that both
   states of a homologous pair fall into the same group. (FFT anchoring
   uses amino acid polarity and volume and has no 3Di counterpart, so
   tools/mafft_3di.sh turns it off.)

The relabeling is stored in the "# 3Di letters:" and "# MAFFT letters:"
comment lines of the output, from which tools/mafft_3di.sh translates
sequences and alignments. Usage:

    python tools/make_mafft_3di_matrix.py tools/mat3di.out \\
        tools/mat3di_mafft.out
"""

import argparse
import math
import random
from pathlib import Path

# MAFFT's amino acid order and its 6-mer groups (locgrpd in core/blosum.c).
MAFFT_ORDER = "ARNDCQEGHILKMFPSTWYV"
MAFFT_GROUPS = ["AGPST", "ILMV", "NDQE", "RHK", "FWY", "C"]


def read_foldseek_matrix(path: Path):
    """Returns letters, scores, background frequencies and lambda."""
    background, lam, rows = None, None, []
    for line in path.read_text().splitlines():
        if line.startswith("# Background"):
            background = [float(x) for x in line.split(":")[1].split()]
        elif line.startswith("# Lambda"):
            lam = float(line.split(":")[1])
        elif line.strip() and not line.startswith("#"):
            rows.append(line.split())
    if background is None or lam is None or not rows:
        raise ValueError(f"{path} is not a foldseek substitution matrix.")
    # Drop the X row and column, which MAFFT does not read.
    letters = [c for c in rows[0] if c != "X"]
    scores = {
        (row[0], col): float(value)
        for row in rows[1:] if row[0] != "X"
        for col, value in zip(rows[0], row[1:]) if col != "X"
    }
    total = sum(background[:len(letters)])
    freq = {c: f / total for c, f in zip(letters, background)}
    return letters, scores, freq, lam


def same_group_mass(groups, pair_prob) -> float:
    """Probability that a homologous pair of distinct states shares a
    group."""
    return sum(
        pair_prob[a, b] for group in groups for a in group for b in group
        if a != b
    )


def partition(letters, pair_prob, restarts: int, seed: int):
    """Partitions the states into groups of the sizes of MAFFT_GROUPS by
    best-improvement swaps from random starts."""
    rng = random.Random(seed)
    sizes = [len(g) for g in MAFFT_GROUPS]
    best, best_mass = None, -1.0
    for _ in range(restarts):
        states = list(letters)
        rng.shuffle(states)
        groups, start = [], 0
        for size in sizes:
            groups.append(states[start:start + size])
            start += size
        mass = same_group_mass(groups, pair_prob)
        improved = True
        while improved:
            improved = False
            for gi in range(len(groups)):
                for gj in range(gi + 1, len(groups)):
                    for i in range(len(groups[gi])):
                        for j in range(len(groups[gj])):
                            gi_state = groups[gi][i]
                            gj_state = groups[gj][j]
                            groups[gi][i], groups[gj][j] = gj_state, gi_state
                            new = same_group_mass(groups, pair_prob)
                            if new > mass + 1e-12:
                                mass, improved = new, True
                            else:
                                groups[gi][i] = gi_state
                                groups[gj][j] = gj_state
        if mass > best_mass:
            best, best_mass = [list(g) for g in groups], mass
    return best, best_mass


def relabel(groups) -> dict[str, str]:
    """Maps each 3Di state to a letter of the MAFFT group it was put in.
    Within a group, states keep their letter where possible."""
    mapping = {}
    for states, letters in zip(groups, MAFFT_GROUPS):
        free = [c for c in letters if c not in states]
        for state in sorted(states):
            mapping[state] = state if state in letters else free.pop(0)
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0]
    )
    parser.add_argument(
        "foldseek_matrix", type=Path,
        help="foldseek's mat3di.out (data/mat3di.out in its repository)."
    )
    parser.add_argument("output", type=Path, help="MAFFT matrix to write.")
    parser.add_argument("--restarts", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    letters, scores, freq, lam = read_foldseek_matrix(args.foldseek_matrix)
    if sorted(letters) != sorted(MAFFT_ORDER):
        parser.error("Expected 20 3Di states named like amino acids.")
    # Implied pair probabilities of the log-odds scores
    # s = log(p_ab / (f_a f_b)) / lambda.
    pair_prob = {
        (a, b): freq[a] * freq[b] * math.exp(lam * scores[a, b])
        for a in letters for b in letters
    }
    norm = sum(pair_prob.values())
    pair_prob = {k: v / norm for k, v in pair_prob.items()}

    groups, mass = partition(letters, pair_prob, args.restarts, args.seed)
    mapping = relabel(groups)
    inverse = {v: k for k, v in mapping.items()}
    diag = sum(pair_prob[a, a] for a in letters)
    print("Groups (3Di states):",
          " ".join("".join(sorted(g)) for g in groups))
    print(f"P(same state) = {diag:.3f}, "
          f"P(distinct states, same group) = {mass:.3f}")

    states = "".join(sorted(letters))
    lines = [
        f"# foldseek 3Di matrix ({args.foldseek_matrix.name}) for MAFFT "
        "--aamatrix, written by tools/make_mafft_3di_matrix.py.",
        "# Rows and columns are named by the MAFFT letters below. Use "
        "tools/mafft_3di.sh, which translates sequences accordingly.",
        f"# 3Di letters: {states}",
        f"# MAFFT letters: {''.join(mapping[c] for c in states)}",
        "   " + "".join(f"{c:>4}" for c in MAFFT_ORDER),
    ]
    for row in MAFFT_ORDER:
        lines.append(f"{row}  " + "".join(
            f"{scores[inverse[row], inverse[col]]:4.0f}"
            for col in MAFFT_ORDER
        ))
    lines.append("frequency " + " ".join(
        f"{freq[inverse[c]]:.6f}" for c in MAFFT_ORDER
    ))
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
