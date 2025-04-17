## Projects a multiple sequence alignment (fasta file) with respect to the sequences
## in the reference file (fasta file).
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

argparser = argparse.ArgumentParser(
    description="Projects a multiple sequence alignment (fasta file)" "" \
    "with respect to the sequences in the reference file (fasta file)."
)

argparser.add_argument(
    "--msa",
    type=str,
    required=True,
    help="Path to the multiple sequence alignment file (fasta format).",
)

argparser.add_argument(
    "--ref",
    type=str,
    required=True,
    help="Path to the reference file (fasta format).",
)

argparser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Path to the output file (fasta format).",
)

args = argparser.parse_args()


# Read the reference sequences
ref_seqs = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))

# Read the MSA sequences
msa_seqs = SeqIO.to_dict(SeqIO.parse(args.msa, "fasta"))

#convert to numpy arrays
ref_seqs_list = []
ref_ids = []
for seq, seq_record in ref_seqs.items():
    ref_seqs_list.append(seq_record.seq)
    ref_ids.append(seq)

msa_seqs_list = []
msa_ids = []
for seq, seq_record in msa_seqs.items():
    msa_seqs_list.append(seq_record.seq)
    msa_ids.append(seq)

ref_seqs_array = np.array(ref_seqs_list)
msa_seqs_array = np.array(msa_seqs_list)

#find the subset
ref_indices = [msa_ids.index(ref_id) for ref_id in ref_ids if ref_id in msa_ids]

projected_seqs_array = msa_seqs_array[ref_indices]

# remove gap only columns
gap_only = np.all((projected_seqs_array == "-") | (projected_seqs_array == "."), axis=0)
projected_seqs_array = projected_seqs_array[:, ~gap_only]

# Write the projected sequences to a new fasta file
with open(args.out, "w") as out_file:
    for i, seq_id in enumerate(ref_ids):
        seq_record = SeqIO.SeqRecord(
            Seq("".join(projected_seqs_array[i]).upper().replace(".", "-")),
            id=seq_id,
            description="",
        )
        SeqIO.write(seq_record, out_file, "fasta")




