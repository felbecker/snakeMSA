import argparse
import gzip

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

GZIP_MAGIC = b"\x1f\x8b"


def open_msa(path):
    """Open an MSA as text, decompressing it on the fly if it is gzipped
    (e.g. learnMSA --compress). Detected by content, not by file name, as
    the pipeline keeps the alignment's name without a ".gz" suffix.
    """
    with open(path, "rb") as f:
        is_gzip = f.read(2) == GZIP_MAGIC
    return gzip.open(path, "rt") if is_gzip else open(path)


def main() -> None:
    parser = argparse.ArgumentParser( 
        description="Project MSA based on reference IDs."
    )
    parser.add_argument(
        "--msa", required=True, help="Path to MSA FASTA."
    )
    parser.add_argument(
        "--ref", required=True, help="Path to reference FASTA."
    )
    parser.add_argument(
        "--out", required=True, help="Output FASTA path."
    )
    args = parser.parse_args()

    # Step 1: Read reference IDs in order
    ref_ids = []
    with open(args.ref) as ref_file:
        for record in SeqIO.parse(ref_file, "fasta"):
            ref_ids.append(record.id)

    # Step 2: Load matching MSA sequences
    ref_seqs = {}
    with open_msa(args.msa) as msa_file:
        for record in SeqIO.parse(msa_file, "fasta"):
            if record.id in ref_ids:
                ref_seqs[record.id] = str(record.seq).upper()
                if len(ref_seqs) == len(ref_ids):  # Early stop if all found
                    break

    # Ensure all references found
    missing = [rid for rid in ref_ids if rid not in ref_seqs]
    if missing:
        raise ValueError(f"Missing reference IDs in MSA: {missing}")

    # Step 3: Convert to NumPy array
    aligned = np.array(
        [list(ref_seqs[rid]) for rid in ref_ids], dtype='U1'
    )  # shape: (n_refs, L)

    # Step 4: Remove gap-only columns
    gaps = (aligned == "-") | (aligned == ".")
    keep_columns = ~np.all(gaps, axis=0)
    projected = aligned[:, keep_columns]

    # Step 5: Write result
    with open(args.out, "w") as out_file:
        for i, rid in enumerate(ref_ids):
            seq_str = "".join(projected[i]).replace(".", "-")
            record = SeqRecord(Seq(seq_str), id=rid, description="")
            SeqIO.write(record, out_file, "fasta")

if __name__ == "__main__":
    main()
