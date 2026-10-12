"""Predict per-residue ProstT5 3Di logits for all FASTA files of a directory.

Runs ``learnMSA.structure.predict_3di`` once per file and stores
``<family>.3di.npz`` (and optionally ``<family>.3di.fasta``) in the given
output directories, e.g. to build data/homfam/prostt5_outputs:

    python tools/predict_prostt5.py data/homfam/unaligned \\
        data/homfam/prostt5_outputs

Run it in an environment with the learnMSA PyTorch extras (learnMSAtorch).
Existing outputs are kept, so an interrupted run can simply be restarted.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

#: Suffixes dropped from the file name to get the family name. Other dots are
#: part of the name.
FASTA_SUFFIXES = (".vie", ".fasta", ".fa", ".fas", ".faa", ".fna")

#: learnMSA checkout next to the snakeMSA root.
DEFAULT_LEARNMSA = Path(__file__).resolve().parents[2] / "learnMSA"


def family_name(path: Path) -> str:
    """File name without a trailing FASTA suffix (none for e.g. 'kinase')."""
    if path.suffix in FASTA_SUFFIXES:
        return path.stem
    return path.name


def collect_inputs(input_dir: Path) -> dict[str, Path]:
    """Maps family name to file for all regular, non-hidden files."""
    inputs: dict[str, Path] = {}
    for path in sorted(input_dir.iterdir()):
        if path.name.startswith(".") or not path.is_file():
            continue
        name = family_name(path)
        if name in inputs:
            raise ValueError(
                f"{inputs[name].name} and {path.name} both map to the family "
                f"name '{name}'."
            )
        inputs[name] = path
    return inputs


def predict(
    fasta: Path,
    npz: Path,
    fasta_out: Path | None,
    learnmsa: Path,
    extra_args: list[str],
) -> None:
    """Runs predict_3di on one file. Outputs appear only on success."""
    # Write under temporary names so that an interrupted run never leaves a
    # truncated output behind that later runs would skip.
    tmp_npz = npz.with_suffix(".tmp.npz")
    tmp_fasta = fasta_out.with_suffix(".tmp.fasta") if fasta_out else None
    cmd = [
        sys.executable, "-m", "learnMSA.structure.predict_3di",
        "-i", str(fasta), "-o", str(tmp_npz),
    ]
    if tmp_fasta is not None:
        cmd += ["--fasta", str(tmp_fasta)]
    cmd += extra_args
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(learnmsa)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
    )
    try:
        subprocess.run(cmd, env=env, check=True)
        # os.replace also overwrites existing read-only files.
        os.replace(tmp_npz, npz)
        npz.chmod(0o444)
        if tmp_fasta is not None and fasta_out is not None:
            os.replace(tmp_fasta, fasta_out)
    finally:
        tmp_npz.unlink(missing_ok=True)
        if tmp_fasta is not None:
            tmp_fasta.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="Further arguments (--device, --cache_dir, --split_length, "
        "--batch_memory, --silent) are passed on to predict_3di.",
    )
    parser.add_argument(
        "input_dir", type=Path,
        help="Directory with amino acid FASTA files (any file name)."
    )
    parser.add_argument(
        "npz_dir", type=Path,
        help="Output directory for the <family>.3di.npz logits."
    )
    parser.add_argument(
        "--fasta_dir", type=Path, default=None,
        help="Also write the argmax 3Di sequences as <family>.3di.fasta "
        "to this directory."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Recompute families whose outputs already exist."
    )
    parser.add_argument(
        "--learnmsa", type=Path, default=DEFAULT_LEARNMSA,
        help="learnMSA checkout providing learnMSA.structure.predict_3di "
        "(default: %(default)s)."
    )
    args, extra_args = parser.parse_known_args()

    inputs = collect_inputs(args.input_dir)
    if not inputs:
        parser.error(f"No files found in {args.input_dir}.")
    args.npz_dir.mkdir(parents=True, exist_ok=True)
    if args.fasta_dir is not None:
        args.fasta_dir.mkdir(parents=True, exist_ok=True)

    done, skipped, failed = [], [], []
    for i, (name, path) in enumerate(inputs.items(), start=1):
        npz = args.npz_dir / f"{name}.3di.npz"
        fasta_out = (
            args.fasta_dir / f"{name}.3di.fasta"
            if args.fasta_dir is not None else None
        )
        outputs = [npz] + ([fasta_out] if fasta_out is not None else [])
        if not args.overwrite and all(o.exists() for o in outputs):
            print(f"[{i}/{len(inputs)}] {name}: skipped (exists)", flush=True)
            skipped.append(name)
            continue
        print(f"[{i}/{len(inputs)}] {name}", flush=True)
        try:
            predict(path, npz, fasta_out, args.learnmsa, extra_args)
            done.append(name)
        except subprocess.CalledProcessError as e:
            print(f"{name}: failed ({e})", file=sys.stderr, flush=True)
            failed.append(name)

    print(
        f"Done: {len(done)}, skipped: {len(skipped)}, failed: {len(failed)}."
    )
    if failed:
        print("Failed: " + " ".join(failed), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
