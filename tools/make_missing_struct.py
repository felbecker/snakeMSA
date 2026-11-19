"""Build a benchmark that simulates missing structures from afdb-clusters-test.

For every cluster and every fraction p, the sample <rep>_p<100p> has all
sequences of the cluster but structures for only round(p * n) of its n
members. The structures kept are nested across p: each cluster gets one
random member order (fixed seed) and every p keeps a prefix of it. Layout of
the output (relative symlinks into the source data, which is read-only):

    unaligned/<sample>.fa   all sequences of the cluster
    pdb/<sample>/           the kept structures (empty for p = 0)
    ref_pdb/<sample>        all structures, so that scoring covers all rows
    structured/<sample>.txt IDs of the kept structures
    manifest.tsv            sample, rep, bin, p, n_members, n_structures

Use it with "pdb": "pdb" and "ref_pdb": "ref_pdb" in the snakeMSA config:

    python tools/make_missing_struct.py
"""

import argparse
import csv
import os
import random
from pathlib import Path

STRUCTURE_SUFFIXES = (".pdb", ".cif", ".ent", ".mmcif")


def read_ids(fasta: Path) -> list[str]:
    with open(fasta) as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def sample_name(rep: str, p: float) -> str:
    # No dots in sample names: snakeMSA splits extensions off them.
    return f"{rep}_p{round(100 * p):02d}"


def link(target: Path, link_path: Path) -> None:
    """A relative symlink to the resolved target."""
    target = target.resolve()
    os.symlink(os.path.relpath(target, link_path.parent), link_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--src", type=Path,
                        default=Path("data/afdb-clusters-test"))
    parser.add_argument("--out", type=Path,
                        default=Path("data/afdb-clusters-missing-struct"))
    parser.add_argument("--p", type=float, nargs="+",
                        default=[0, 0.1, 0.25, 0.5, 0.75])
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.out.exists():
        raise SystemExit(f"{args.out} exists; remove it to rebuild.")
    if any(not 0 <= p <= 1 for p in args.p):
        raise SystemExit("--p values must be in [0, 1].")
    with open(args.src / "split.tsv") as f:
        bins = {r["rep_accession"]: r["bin"]
                for r in csv.DictReader(f, delimiter="\t")}

    subdirs = ("unaligned", "pdb", "ref_pdb", "structured")
    for sub in subdirs:
        (args.out / sub).mkdir(parents=True)
    rows = []
    for fasta in sorted((args.src / "unaligned").glob("*.fa")):
        rep = fasta.stem
        pdb_dir = args.src / "pdb" / rep
        files = {
            path.stem: path.name for path in pdb_dir.iterdir()
            if path.suffix.lower() in STRUCTURE_SUFFIXES
        }
        ids = read_ids(fasta)
        missing = sorted(set(ids) - files.keys())
        if missing:
            raise SystemExit(
                f"{rep}: {len(missing)} members without a structure, "
                f"e.g. {missing[0]}."
            )
        order = sorted(ids)
        random.Random(f"{args.seed}:{rep}").shuffle(order)
        for p in args.p:
            sample = sample_name(rep, p)
            k = int(p * len(order) + 0.5)
            link(fasta, args.out / "unaligned" / f"{sample}.fa")
            (args.out / "pdb" / sample).mkdir()
            for seq_id in order[:k]:
                name = files[seq_id]
                link(pdb_dir / name, args.out / "pdb" / sample / name)
            link(pdb_dir, args.out / "ref_pdb" / sample)
            with open(args.out / "structured" / f"{sample}.txt", "w") as f:
                f.writelines(f"{seq_id}\n" for seq_id in order[:k])
            rows.append((sample, rep, bins.get(rep, ""), p, len(order), k))

    with open(args.out / "manifest.tsv", "w") as f:
        f.write("sample\trep\tbin\tp\tn_members\tn_structures\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"{args.out}: {len(rows)} samples from "
          f"{len(rows) // len(args.p)} clusters, p = {args.p}")


if __name__ == "__main__":
    main()
