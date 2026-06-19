"""Parses a tsv files with header "msa_name", "lddt_foldmason",
"lddt_foldmason_core", "lddt_muscle" and computes the average LDDT scores
across all MSAs."""

import argparse
import csv

argparser = argparse.ArgumentParser(
    description="Parses a tsv files with header 'msa_name', 'lddt_foldmason', " \
    "'lddt_foldmason_core', 'lddt_muscle' and computes the average LDDT scores across all MSAs."
)
argparser.add_argument(
    "--input", type=str,
    default="lddt_scores.tsv",
    help="Path to the input TSV file."
)
args = argparser.parse_args()

if __name__ == '__main__':

    input_path = args.input

    total_foldmason = 0.0
    total_foldmason_core = 0.0
    total_muscle = 0.0
    count = 0

    with open(input_path, "r", newline="") as in_f:
        tsv_reader = csv.DictReader(in_f, delimiter="\t")
        for row in tsv_reader:
            total_foldmason += float(row["lddt_foldmason"])
            total_foldmason_core += float(row["lddt_foldmason_core"])
            total_muscle += float(row["lddt_muscle"])
            count += 1

    print(f"Average LDDT FoldMason: {total_foldmason / count:.4f}")
    print(f"Average LDDT FoldMason Core: {total_foldmason_core / count:.4f}")
    print(f"Average LDDT Muscle: {total_muscle / count:.4f}")
