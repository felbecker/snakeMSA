"""Predict per-residue ProstT5 3Di logits for all FASTA files of a directory.

Runs ``learnMSA.structure.predict_3di`` once per file and stores
``<family>.3di.npz`` (and optionally ``<family>.3di.fasta``) in the given
output directories, e.g. to build data/homfam/prostt5_outputs:

    python tools/predict_prostt5.py data/homfam/unaligned \\
        data/homfam/prostt5_outputs

With several visible GPUs, one predict job per GPU runs in parallel (each job
sees only its own GPU), so the output of the jobs may interleave.

Run it in an environment with the learnMSA PyTorch extras (learnMSAtorch).
Existing outputs are kept, so an interrupted run can simply be restarted.
"""

import argparse
import os
import queue
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
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


def parse_gpus(spec: str) -> list[str]:
    """Splits a comma-separated GPU list; '' and '-1' mean no GPUs."""
    gpus = [g.strip() for g in spec.split(",") if g.strip()]
    return [] if gpus == ["-1"] else gpus


def visible_gpus() -> list[str]:
    """GPU ids from CUDA_VISIBLE_DEVICES or, if unset, from nvidia-smi."""
    if "CUDA_VISIBLE_DEVICES" in os.environ:
        return parse_gpus(os.environ["CUDA_VISIBLE_DEVICES"])
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
            capture_output=True, text=True, check=True,
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError):
        return []
    return [line.strip() for line in out.splitlines() if line.strip()]


def predict(
    fasta: Path,
    npz: Path,
    fasta_out: Path | None,
    learnmsa: Path,
    extra_args: list[str],
    gpu: str | None = None,
) -> None:
    """Runs predict_3di on one file. Outputs appear only on success.

    If ``gpu`` is given, the job sees only that GPU.
    """
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
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = gpu
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
        "--batch_memory, --silent) are passed on to predict_3di. With "
        "--device, a single job runs on that device.",
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
    parser.add_argument(
        "--gpus", type=parse_gpus, default=None,
        help="Comma-separated GPU ids, one parallel job per GPU (default: "
        "CUDA_VISIBLE_DEVICES, else all GPUs listed by nvidia-smi)."
    )
    args, extra_args = parser.parse_known_args()

    inputs = collect_inputs(args.input_dir)
    if not inputs:
        parser.error(f"No files found in {args.input_dir}.")
    args.npz_dir.mkdir(parents=True, exist_ok=True)
    if args.fasta_dir is not None:
        args.fasta_dir.mkdir(parents=True, exist_ok=True)

    gpus = visible_gpus() if args.gpus is None else args.gpus
    if any(a == "--device" or a.startswith("--device=") for a in extra_args):
        if len(gpus) > 1:
            print("--device given: running a single job.", flush=True)
        gpus = []
    print(
        f"Running {max(1, len(gpus))} job(s) in parallel"
        + (f" on GPU(s) {','.join(gpus)}." if gpus else "."),
        flush=True,
    )

    done, skipped, failed = [], [], []
    todo = []
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
        todo.append((i, name, path, npz, fasta_out))

    # Hands each running job one free GPU (None without GPUs).
    free_gpus: queue.Queue[str | None] = queue.Queue()
    for gpu in gpus or [None]:
        free_gpus.put(gpu)

    def run(i, name, path, npz, fasta_out):
        gpu = free_gpus.get()
        try:
            on_gpu = f" (GPU {gpu})" if gpu is not None else ""
            print(f"[{i}/{len(inputs)}] {name}{on_gpu}", flush=True)
            predict(path, npz, fasta_out, args.learnmsa, extra_args, gpu)
        finally:
            free_gpus.put(gpu)

    executor = ThreadPoolExecutor(max_workers=max(1, len(gpus)))
    try:
        futures = {executor.submit(run, *job): job[1] for job in todo}
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
                done.append(name)
            except subprocess.CalledProcessError as e:
                print(f"{name}: failed ({e})", file=sys.stderr, flush=True)
                failed.append(name)
    except KeyboardInterrupt:
        # The jobs share our process group and receive the SIGINT as well.
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    executor.shutdown()

    print(
        f"Done: {len(done)}, skipped: {len(skipped)}, failed: {len(failed)}."
    )
    if failed:
        print("Failed: " + " ".join(failed), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
