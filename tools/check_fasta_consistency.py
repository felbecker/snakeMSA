"""Checks if two FASTA files have the same sequence IDs and lengths."""

# Example usage:
# for file in data/homfam/unaligned/*.vie; do echo $file && python tools/check_fasta_consistency.py $file data/homfam/3di/$(basename "$file" .vie).fasta; done

import argparse

from learnMSA.util import SequenceDataset


argparser = argparse.ArgumentParser(
    description="Checks if two FASTA files have the same sequence IDs and " \
    "lengths."
)
argparser.add_argument("fasta_1", type=str, help="First FASTA file.")
argparser.add_argument("fasta_2", type=str, help="Second FASTA file.")
argparser.add_argument(
    "--ordered", action="store_true",
    help="Also check if sequences are in the same order."
)

args = argparser.parse_args()

if __name__ == '__main__':

    fasta_1 = SequenceDataset(args.fasta_1)
    fasta_2 = SequenceDataset(args.fasta_2)

    assert len(fasta_1) == len(fasta_2), \
        "FASTA files have different number of sequences."
    ids_1, ids_2 = set(fasta_1.seq_ids), set(fasta_2.seq_ids)
    only_in_1 = ids_1 - ids_2
    only_in_2 = ids_2 - ids_1
    assert not only_in_1 and not only_in_2, (
        f"FASTA files have different sequence IDs. "
        f"Only in fasta_1: {only_in_1}. Only in fasta_2: {only_in_2}."
    )
    mismatched = [
        fasta_1.seq_ids[i]
        for i in range(len(fasta_1))
        if fasta_1.seq_lens[i] != fasta_2.seq_lens[fasta_2.index(fasta_1.seq_ids[i])]
    ]
    assert not mismatched, (
        f"FASTA files have different sequence lengths for: {mismatched}."
    )

    if args.ordered:
        assert fasta_1.seq_ids == fasta_2.seq_ids, \
            "FASTA files have different sequence order."
        assert all(fasta_1.seq_lens == fasta_2.seq_lens), \
            "FASTA files have different sequence lengths in the same order."

    print("FASTA files are consistent.")