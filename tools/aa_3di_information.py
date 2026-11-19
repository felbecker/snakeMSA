#!/usr/bin/env python3
"""Information content of amino acids (aa) and 3Di at alignment columns.

Two analyses that back the choice of conditionally independent aa and 3Di
emissions in learnMSA. All numbers are in bits.

``hmm``
    Entropy and information content (KL divergence to the background) of the
    aa and 3Di match-state emission distributions of trained learnMSA models
    (selected head only). Needs torch to read the ``.model.pt`` checkpoints,
    e.g. ``conda run -n learnMSA python tools/aa_3di_information.py hmm``.

``homstrad``
    Mutual information between aa (A), 3Di (T) and the column (C) of
    structure-based HOMSTRAD alignments with 3Di computed from the PDB
    structures. The information content of a column under a joint emission
    model exceeds that of an independent model by exactly I(A;T|C):

        E_c KL(p(a,t|c) || q(a)q(t)) = I(A;C) + I(T;C) + I(A;T|C),

    where q are the family's marginal frequencies. I(A;T|C) is also the
    expected per-residue log-likelihood gain of a perfectly fitted joint
    emitter over a perfectly fitted independent one. Plug-in estimates are
    corrected by subtracting the mean of permutation nulls. Numpy only.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

AA_ALPHABET = "ARNDCQEGHILKMFPSTWYV"
TDI_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
GAPS = "-."
D = 20


# ---------------------------------------------------------------------------
# Information-theoretic helpers (bits)
# ---------------------------------------------------------------------------

def entropy_rows(p: np.ndarray) -> np.ndarray:
    """Entropy (bits) of each row of a probability matrix."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lp = np.where(p > 0, np.log2(p), 0.0)
    return -(p * lp).sum(-1)


def kl_rows(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """KL(p_i || q) in bits for each row p_i."""
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.where(p > 0, np.log2(p / q), 0.0)
    return (p * lr).sum(-1)


def entropy_counts(counts: np.ndarray) -> float:
    """Plug-in entropy (bits) of a count vector of any shape."""
    c = counts[counts > 0].astype(np.float64)
    n = c.sum()
    return float(-(c / n * np.log2(c / n)).sum())


def cond_entropy(group: np.ndarray, x: np.ndarray, k: int) -> float:
    """Plug-in H(X | G) in bits; ``x`` takes values in [0, k)."""
    n = len(x)
    gc = np.bincount(group)
    joint = np.bincount(group * k + x, minlength=len(gc) * k)
    joint = joint[joint > 0].astype(np.float64)
    gc = gc[gc > 0].astype(np.float64)
    return float(
        (-(joint * np.log2(joint)).sum() + (gc * np.log2(gc)).sum()) / n
    )


def mi_given(group: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    """Plug-in I(X;Y|G) = H(X|G) + H(Y|G) - H(X,Y|G) in bits."""
    return (cond_entropy(group, x, D) + cond_entropy(group, y, D)
            - cond_entropy(group, x * D + y, D * D))


def mi_with_group(group: np.ndarray, x: np.ndarray) -> float:
    """Plug-in I(X;G) = H(X) - H(X|G) in bits."""
    return entropy_counts(np.bincount(x, minlength=D)) \
        - cond_entropy(group, x, D)


def mi_table(table: np.ndarray) -> float:
    """Plug-in I(X;Y) in bits from a 2D count table."""
    return (entropy_counts(table.sum(1)) + entropy_counts(table.sum(0))
            - entropy_counts(table))


def shuffle_within(group: np.ndarray, x: np.ndarray,
                   rng: np.random.Generator) -> np.ndarray:
    """Permute ``x`` within each group. ``group`` must be sorted."""
    order = np.lexsort((rng.random(len(x)), group))
    return x[order]


# ---------------------------------------------------------------------------
# hmm: entropy of trained match states
# ---------------------------------------------------------------------------

def softmax_rows(x: np.ndarray) -> np.ndarray:
    x = x - x.max(-1, keepdims=True)
    e = np.exp(x)
    return e / e.sum(-1, keepdims=True)


def head_matrices(kernel: np.ndarray, lengths) -> list[np.ndarray]:
    """Split a flat learnMSA profile-emitter kernel into per-head matrices.

    Per head the kernel holds (L+1) x 20 logits: L match rows followed by
    one insert row. Checked against ``TorchProfileEmitter.matrix()``.
    """
    mats, off = [], 0
    for L in lengths:
        n = (L + 1) * D
        mats.append(softmax_rows(kernel[off:off + n].reshape(L + 1, D)))
        off += n
    assert off == len(kernel), "kernel size does not match model lengths"
    return mats


def selected_head(log_file: str) -> int:
    best = None
    with open(log_file) as f:
        for line in f:
            m = re.search(r"Best model:\s+(\d+)", line)
            if m:
                best = int(m.group(1))
    if best is None:
        raise ValueError(f"no 'Best model' line in {log_file}")
    return best


def run_hmm(args: argparse.Namespace) -> None:
    import torch  # only needed here

    rows = []
    for run_dir in sorted(glob.glob(os.path.join(args.outputs, args.runs))):
        run = os.path.basename(run_dir)
        for wd in sorted(glob.glob(os.path.join(run_dir, "workdirs", "*"))):
            name = os.path.basename(wd)
            fam = os.path.splitext(name)[0]
            ckpt = os.path.join(wd, f"{fam}.model.pt")
            log = os.path.join(run_dir, "logs", f"{name}.log")
            if not os.path.exists(ckpt) or not os.path.exists(log):
                print(f"skip {run}/{name}: missing model or log",
                      file=sys.stderr)
                continue
            d = torch.load(ckpt, map_location="cpu", weights_only=False)
            lengths = list(d["context"]["model_lengths"])
            sd = d["state_dict"]
            h = selected_head(log)
            stats = {}
            for track, key in (("aa", "profile_emitter"),
                               ("3di", "struct_emitter")):
                k = sd[f"phmm_layer.{key}.kernel"].double().numpy()
                M = head_matrices(k, lengths)[h]
                L = lengths[h]
                match, background = M[:L], M[L]
                stats[f"H_{track}"] = entropy_rows(match)
                stats[f"IC_{track}"] = kl_rows(match, background)
                stats[f"Hbg_{track}"] = float(entropy_rows(background))
            rows.append(dict(
                run=run, family=fam, L=lengths[h],
                H_aa=stats["H_aa"].mean(), H_3di=stats["H_3di"].mean(),
                IC_aa=stats["IC_aa"].mean(), IC_3di=stats["IC_3di"].mean(),
                frac_H3di_lt_Haa=float(
                    (stats["H_3di"] < stats["H_aa"]).mean()),
                frac_IC3di_gt_ICaa=float(
                    (stats["IC_3di"] > stats["IC_aa"]).mean()),
                Hbg_aa=stats["Hbg_aa"], Hbg_3di=stats["Hbg_3di"],
            ))

    if not rows:
        sys.exit("no models found")
    cols = list(rows[0].keys())
    with open(args.tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    runs = sorted({r["run"] for r in rows})
    metrics = ["H_aa", "H_3di", "IC_aa", "IC_3di",
               "frac_H3di_lt_Haa", "frac_IC3di_gt_ICaa", "Hbg_aa", "Hbg_3di"]
    print(f"# hmm: {len(rows)} models, {len(runs)} runs, "
          f"{len({r['family'] for r in rows})} families "
          "(selected head, match states; bits)")
    print("# per run: mean over families (median in brackets)")
    per_run = {m: [] for m in metrics}
    for run in runs:
        rr = [r for r in rows if r["run"] == run]
        line = [run, f"n={len(rr)}"]
        for m in metrics:
            v = np.array([r[m] for r in rr])
            per_run[m].append(v.mean())
            line.append(f"{m}={v.mean():.3f} [{np.median(v):.3f}]")
        print("  ".join(line))
    print("# across runs: mean +- SD of the per-run family means")
    for m in metrics:
        v = np.array(per_run[m])
        print(f"{m}\t{v.mean():.3f}\t+- {v.std(ddof=1) if len(v) > 1 else 0:.3f}")
    fam_ok = {}
    for r in rows:
        fam_ok.setdefault(r["family"], []).append(r["H_3di"] < r["H_aa"])
    n_lower = sum(all(v) for v in fam_ok.values())
    print(f"# families with lower mean 3Di than aa entropy in all runs: "
          f"{n_lower}/{len(fam_ok)}")
    fam_ic = {}
    for r in rows:
        fam_ic.setdefault(r["family"], []).append(r["IC_3di"] > r["IC_aa"])
    n_ic = sum(all(v) for v in fam_ic.values())
    print(f"# families with higher mean 3Di than aa IC in all runs: "
          f"{n_ic}/{len(fam_ic)}")
    print(f"# per-family table written to {args.tsv}")


# ---------------------------------------------------------------------------
# homstrad: I(A;C), I(T;C), I(A;T|C), I(A;T) on structural alignments
# ---------------------------------------------------------------------------

def read_fasta(path: str) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def encode(s: str, alphabet: str) -> np.ndarray:
    """Letters -> [0, 20); gaps -> -1; unknown letters -> -2."""
    lut = np.full(256, -2, dtype=np.int64)
    for i, c in enumerate(alphabet):
        lut[ord(c)] = i
        lut[ord(c.lower())] = i
    for g in GAPS:
        lut[ord(g)] = -1
    return lut[np.frombuffer(s.encode(), dtype=np.uint8)]


def load_family(aa_file: str, tdi_file: str, min_occ: float):
    """Return (column, aa, 3di) arrays for residues in match-like columns.

    Residues are sorted by column. Columns with a non-gap fraction below
    ``min_occ`` are dropped, as are residues with an unknown letter.
    """
    aa, tdi = read_fasta(aa_file), read_fasta(tdi_file)
    if set(aa) != set(tdi):
        raise ValueError(f"{aa_file}: sequence IDs differ from {tdi_file}")
    ids = sorted(aa)
    A = np.stack([encode(aa[i], AA_ALPHABET) for i in ids])
    T = np.stack([encode(tdi[i], TDI_ALPHABET) for i in ids])
    if A.shape != T.shape or not np.array_equal(A == -1, T == -1):
        raise ValueError(f"{aa_file}: gap patterns of aa and 3Di differ")
    occ = (A != -1).mean(0)
    keep_cols = np.where(occ >= min_occ)[0]
    A, T = A[:, keep_cols], T[:, keep_cols]
    col = np.broadcast_to(np.arange(A.shape[1]), A.shape)
    ok = (A >= 0) & (T >= 0)
    # Flatten column-major so that residues are grouped by column.
    col, a, t = col.T[ok.T], A.T[ok.T], T.T[ok.T]
    return col.astype(np.int64), a, t, len(ids), len(keep_cols)


def family_stats(col, a, t, B, rng):
    """Plug-in and permutation-null mean for each quantity of one family."""
    plug = np.array([mi_with_group(col, a), mi_with_group(col, t),
                     mi_given(col, a, t)])
    null = np.zeros(3)
    for _ in range(B):
        null += np.array([
            mi_with_group(col, rng.permutation(a)),
            mi_with_group(col, rng.permutation(t)),
            mi_given(col, a, shuffle_within(col, t, rng)),
        ])
    return plug, null / B


def run_homstrad(args: argparse.Namespace) -> None:
    rng = np.random.default_rng(args.seed)
    fams, per_fam = [], []
    table = np.zeros((D, D), dtype=np.int64)
    n_cols_total, col_sizes = 0, []
    for tdi_file in sorted(glob.glob(os.path.join(args.tdi_dir, "*.fasta"))):
        fam = os.path.splitext(os.path.basename(tdi_file))[0]
        aa_file = os.path.join(args.aa_dir, os.path.basename(tdi_file))
        if not os.path.exists(aa_file):
            print(f"skip {fam}: no aa alignment", file=sys.stderr)
            continue
        col, a, t, nseq, ncol = load_family(aa_file, tdi_file, args.min_occ)
        if args.shuffle_control:
            t = shuffle_within(col, t, rng)
        plug, null = family_stats(col, a, t, args.permutations, rng)
        tab = np.zeros((D, D), dtype=np.int64)
        np.add.at(tab, (a, t), 1)
        table += tab
        n_cols_total += ncol
        col_sizes.append(np.bincount(col))
        fams.append(fam)
        per_fam.append(dict(family=fam, nseq=nseq, ncol=ncol, n=len(a),
                            plug=plug, null=null, corr=plug - null))
    if not per_fam:
        sys.exit("no families found")

    names = ["I(A;C)", "I(T;C)", "I(A;T|C)"]
    w = np.array([f["n"] for f in per_fam], dtype=np.float64)
    corr = np.stack([f["corr"] for f in per_fam])
    plug = np.stack([f["plug"] for f in per_fam])
    null = np.stack([f["null"] for f in per_fam])

    # Global aa-3Di coupling, pooled over all families.
    i_at = mi_table(table)
    a_all = np.repeat(np.arange(D * D) // D, table.ravel())
    t_all = np.repeat(np.arange(D * D) % D, table.ravel())
    null_at = np.mean([
        mi_table(np.bincount(a_all * D + rng.permutation(t_all),
                             minlength=D * D).reshape(D, D))
        for _ in range(args.permutations)])
    i_at_corr = i_at - null_at

    def summary(idx):
        ww = w[idx]
        v = (corr[idx] * ww[:, None]).sum(0) / ww.sum()
        ic_ind = v[0] + v[1]
        return np.concatenate([v, [ic_ind, ic_ind + v[2], v[2] / ic_ind]])

    point = summary(np.arange(len(per_fam)))
    boot = np.stack([
        summary(rng.integers(0, len(per_fam), len(per_fam)))
        for _ in range(args.bootstrap)])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)

    sizes = np.concatenate(col_sizes)
    print(f"# homstrad: {len(per_fam)} families, {int(w.sum())} residues, "
          f"{n_cols_total} columns with occupancy >= {args.min_occ}; "
          f"residues per column: median {np.median(sizes):.0f}, "
          f"mean {sizes.mean():.1f}")
    if args.shuffle_control:
        print("# CONTROL: 3Di shuffled within columns; I(A;T|C) should be ~0")
    print(f"# residue-weighted means over families, bits; "
          f"{args.permutations} permutations, {args.bootstrap} bootstrap "
          "resamples of families")
    print("quantity\tplug-in\tnull\tcorrected\t95% CI")
    for k, n in enumerate(names):
        pm = (plug[:, k] * w).sum() / w.sum()
        nm = (null[:, k] * w).sum() / w.sum()
        print(f"{n}\t{pm:.4f}\t{nm:.4f}\t{point[k]:.4f}\t"
              f"[{lo[k]:.4f}, {hi[k]:.4f}]")
    print(f"I(A;T) global\t{i_at:.4f}\t{null_at:.4f}\t{i_at_corr:.4f}\t-")
    print(f"IC_indep = I(A;C)+I(T;C)\t\t\t{point[3]:.4f}\t"
          f"[{lo[3]:.4f}, {hi[3]:.4f}]")
    print(f"IC_joint = IC_indep+I(A;T|C)\t\t\t{point[4]:.4f}\t"
          f"[{lo[4]:.4f}, {hi[4]:.4f}]")
    print(f"I(A;T|C) / IC_indep\t\t\t{point[5]:.4f}\t"
          f"[{lo[5]:.4f}, {hi[5]:.4f}]")
    print(f"I(A;T|C) / I(T;C)\t\t\t{point[2] / point[1]:.4f}\t-")
    n_pos = int((corr[:, 2] > 0).sum())
    n_t_gt_a = int((corr[:, 1] > corr[:, 0]).sum())
    print(f"# families with corrected I(A;T|C) > 0: {n_pos}/{len(per_fam)}")
    print(f"# families with I(T;C) > I(A;C): {n_t_gt_a}/{len(per_fam)}")

    with open(args.tsv, "w") as f:
        f.write("family\tnseq\tncol\tresidues\t"
                + "\t".join(f"{n}_plugin\t{n}_corrected" for n in names)
                + "\n")
        for r in per_fam:
            f.write(f"{r['family']}\t{r['nseq']}\t{r['ncol']}\t{r['n']}\t"
                    + "\t".join(f"{p:.4f}\t{c:.4f}"
                                for p, c in zip(r["plug"], r["corr"]))
                    + "\n")
    print(f"# per-family table written to {args.tsv}")


# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    h = sub.add_parser("hmm", help="entropy/IC of trained match states")
    h.add_argument("--outputs",
                   default="outputs/learnMSA_pub_homfam",
                   help="benchmark output directory")
    h.add_argument("--runs", default="learnMSA_struct_replica_*",
                   help="glob of run directories inside --outputs")
    h.add_argument("--tsv", default="hmm_match_entropy.tsv")
    h.set_defaults(func=run_hmm)

    s = sub.add_parser("homstrad", help="MI on HOMSTRAD alignments")
    s.add_argument("--aa-dir", default="data/homstrad/aligned")
    s.add_argument("--tdi-dir", default="data/homstrad/aligned_3Di")
    s.add_argument("--min-occ", type=float, default=0.5,
                   help="minimum non-gap fraction of a column")
    s.add_argument("--permutations", type=int, default=100)
    s.add_argument("--bootstrap", type=int, default=1000)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--shuffle-control", action="store_true",
                   help="shuffle 3Di within columns first (sanity check)")
    s.add_argument("--tsv", default="homstrad_aa_3di_information.tsv")
    s.set_defaults(func=run_homstrad)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
