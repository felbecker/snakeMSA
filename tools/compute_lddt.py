"""Computes different LDDT scores over several alignments and pdb databases."""


import argparse
import csv
import subprocess

from pathlib import Path
import sys


argparser = argparse.ArgumentParser(
    description="Checks if the reference sequences in the homfam collection " \
    "have exact matches to reference sequences in the homstrad collection."
)
argparser.add_argument(
    "--msas", type=str,
    default="data/homfam/aligned",
    help="Path to the msas collection."
)
argparser.add_argument(
    "--pdbs", type=str,
    default="data/homstrad/pdbs",
    help="Path to the pdb files. Should contain a model for each MSA. Sequences" \
    "in the MSA should be identical to the sequences in the pdb files."
)
argparser.add_argument(
    "--workdir", type=str,
    default="tmp",
    help="Path to the working directory."
)
argparser.add_argument(
    "--output", type=str,
    default="lddt_scores.tsv",
    help="Path to the output TSV file."
)
args = argparser.parse_args()



def _parse_foldmason_ouput(output: str) -> str:
    lddt_match = next((
        line for line in output.splitlines()
        if line.startswith("Average MSA LDDT:")
    ), None)
    return lddt_match.split(":")[1].strip() if lddt_match else "-1"


if __name__ == '__main__':

    msa_path = Path(args.msas)
    workdir = Path(args.workdir)
    output_path = Path(args.output)

    workdir.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as out_f:
        tsv_writer = csv.writer(out_f, delimiter="\t")
        tsv_writer.writerow([
            "msa_name", "lddt_foldmason", "lddt_foldmason_core", "lddt_muscle"
        ])

        for msa_file in msa_path.glob("*.vie"):
            msa_filename = msa_file.stem
            print(f"Processing {msa_filename}...")
            pdb_dir = Path(args.pdbs) / msa_filename
            if not pdb_dir.exists():
                print(f"Could not find pdb directory for {msa_file.name}.")
                continue

            # Prepare Foldmason database
            foldmason_db = str(workdir / f"{msa_filename}_db")

            if not Path(foldmason_db).exists():
                cmd = [
                    "foldmason", "createdb",
                    str(pdb_dir),
                    foldmason_db
                ]
                subprocess.run(cmd, check=True)

            try:
                # Foldmason LDDT
                cmd = [
                    "foldmason", "msa2lddt",
                    foldmason_db,
                    str(msa_file)
                ]
                foldmason_out = subprocess.run(
                    cmd, check=True, capture_output=True, text=True
                )
                lddt_foldmason = _parse_foldmason_ouput(foldmason_out.stdout)

                # # Foldmason LDDT core
                cmd = [
                    "foldmason", "msa2lddt",
                    foldmason_db,
                    str(msa_file),
                    "--pair-threshold", "1.0",
                    "--only-scoring-cols", "true"
                ]
                foldmason_out = subprocess.run(
                    cmd, check=True, capture_output=True, text=True
                )
                lddt_foldmason_core = _parse_foldmason_ouput(foldmason_out.stdout)

                # LDDT muscle
                with open(workdir / f"{msa_filename}_pdbfiles.txt", "w") as f:
                    for pdb_file in pdb_dir.glob("*"):
                        f.write(str(pdb_file) + "\n")
                cmd = [
                    "python", "pylddt/src/lddt_mu.py",
                    "--msa", str(msa_file),
                    "--pdbfiles", str(workdir / f"{msa_filename}_pdbfiles.txt"),
                ]
                muscle_out = subprocess.run(
                    cmd, check=True, capture_output=True, text=True
                )
                muscle_lddt = muscle_out.stdout.strip().split("\t")[1]
            except subprocess.CalledProcessError as e:
                print(f"Error processing {msa_filename}: {e}")
                print(f"stdout: {e.stdout}")
                print(f"stderr: {e.stderr}")
                print("Skipping this MSA.")
                continue

            # Write the row for this MSA immediately so partial results
            # are preserved even if a later iteration fails.
            tsv_writer.writerow([
                msa_filename, lddt_foldmason, lddt_foldmason_core, muscle_lddt
            ])
            out_f.flush()

    print(f"Wrote LDDT scores to {output_path}")