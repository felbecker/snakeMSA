"""Translate structures to per-residue 3Di logits for all families of a
directory.

Runs ``learnMSA.structure.translate_3di`` (foldseek's 3Di encoder, no
ProstT5) once per family subdirectory of PDB/mmCIF files and stores
``<family>.3di.npz`` (and optionally ``<family>.3di.fasta`` and the amino
acids of the same chains as ``<family>.fasta``) in the given output
directories, e.g. for data/homstrad/pdb/<family>/*.pdb:

    python tools/translate_3di.py data/homstrad/pdb \\
        data/homstrad/translate_3di_outputs

The .npz has the same format as that of predict_prostt5.py and is marked
as translated from structures, so learnMSA picks the temperature for
structure input (--struct_input_temperature overrides it). Sequence IDs
follow foldseek (file name without extension, plus a chain suffix for
multi-chain files), so use the --aa_fasta_dir output as the matching
amino acid input. Run it in an environment with numpy and gemmi
(e.g. learnMSA). Existing outputs are kept, so an interrupted run can
simply be restarted.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from predict_prostt5 import DEFAULT_LEARNMSA


def collect_families(input_dir: Path) -> dict[str, Path]:
    """Maps family name to directory for all non-hidden subdirectories."""
    return {
        path.name: path
        for path in sorted(input_dir.iterdir())
        if path.is_dir() and not path.name.startswith(".")
    }


def translate(
    structures: Path,
    npz: Path,
    fasta_out: Path | None,
    aa_out: Path | None,
    learnmsa: Path,
    extra_args: list[str],
) -> None:
    """Runs translate_3di on one family. Outputs appear only on success."""
    # Write under temporary names so that an interrupted run never leaves a
    # truncated output behind that later runs would skip.
    outputs = [(npz, npz.with_suffix(".tmp.npz"))]
    cmd = [
        sys.executable, "-m", "learnMSA.structure.translate_3di",
        "-i", str(structures), "-o", str(outputs[0][1]),
    ]
    for flag, path in (("--fasta", fasta_out), ("--aa_fasta", aa_out)):
        if path is not None:
            tmp = path.with_suffix(".tmp.fasta")
            outputs.append((path, tmp))
            cmd += [flag, str(tmp)]
    cmd += extra_args
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(learnmsa)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
    )
    try:
        subprocess.run(cmd, env=env, check=True)
        for final, tmp in outputs:
            # os.replace also overwrites existing read-only files.
            os.replace(tmp, final)
        npz.chmod(0o444)
    finally:
        for _, tmp in outputs:
            tmp.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0].replace("\n", " "),
        epilog="Further arguments (--cache_dir, --silent) are "
        "passed on to translate_3di.",
    )
    parser.add_argument(
        "input_dir", type=Path,
        help="Directory with one subdirectory of structure files per family."
    )
    parser.add_argument(
        "npz_dir", type=Path,
        help="Output directory for the <family>.3di.npz logits."
    )
    parser.add_argument(
        "--fasta_dir", type=Path, default=None,
        help="Also write the 3Di sequences (identical to foldseek) as "
        "<family>.3di.fasta to this directory."
    )
    parser.add_argument(
        "--aa_fasta_dir", type=Path, default=None,
        help="Also write the amino acid sequences of the same chains as "
        "<family>.fasta to this directory."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Recompute families whose outputs already exist."
    )
    parser.add_argument(
        "--learnmsa", type=Path, default=DEFAULT_LEARNMSA,
        help="learnMSA checkout providing learnMSA.structure.translate_3di "
        "(default: %(default)s)."
    )
    args, extra_args = parser.parse_known_args()

    families = collect_families(args.input_dir)
    if not families:
        parser.error(f"No family subdirectories found in {args.input_dir}.")
    for out_dir in (args.npz_dir, args.fasta_dir, args.aa_fasta_dir):
        if out_dir is not None:
            out_dir.mkdir(parents=True, exist_ok=True)

    done, skipped, failed = [], [], []
    for i, (name, path) in enumerate(families.items(), start=1):
        npz = args.npz_dir / f"{name}.3di.npz"
        fasta_out = (
            args.fasta_dir / f"{name}.3di.fasta"
            if args.fasta_dir is not None else None
        )
        aa_out = (
            args.aa_fasta_dir / f"{name}.fasta"
            if args.aa_fasta_dir is not None else None
        )
        outputs = [o for o in (npz, fasta_out, aa_out) if o is not None]
        if not args.overwrite and all(o.exists() for o in outputs):
            print(f"[{i}/{len(families)}] {name}: skipped (exists)",
                  flush=True)
            skipped.append(name)
            continue
        print(f"[{i}/{len(families)}] {name}", flush=True)
        try:
            translate(path, npz, fasta_out, aa_out, args.learnmsa,
                      extra_args)
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
