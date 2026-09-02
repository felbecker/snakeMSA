"""Rewrite reseek/muscle MSA headers back to the input PDB file stems.
"""

import argparse
import os
import re
import sys


def fix(msa_path, pdb_dir, out_path):
    """Map MSA headers onto the PDB file stems in pdb_dir.

    Returns the list of headers that could not be matched (empty = clean).
    """
    stems = {
        os.path.splitext(name)[0]
        for name in os.listdir(pdb_dir)
    }

    unmatched = []
    with open(msa_path) as msa_file, open(out_path, "w") as out_file:
        for line in msa_file:
            if not line.startswith(">"):
                out_file.write(line)
                continue
            fields = line[1:].split()
            token = fields[0] if fields else ""
            if token not in stems:
                stripped = re.sub(r"_[A-Za-z0-9]+$", "", token)
                if stripped in stems:
                    token = stripped
                else:
                    unmatched.append(token)
            out_file.write(">" + token + "\n")
    return unmatched


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rewrite MSA headers to the input PDB file stems."
    )
    parser.add_argument(
        "--msa", required=True, help="Path to MSA FASTA (rewritten in place)."
    )
    parser.add_argument(
        "--pdb", required=True, help="Directory holding the input PDB files."
    )
    parser.add_argument(
        "--out", default=None,
        help="Output FASTA path (default: overwrite --msa)."
    )
    args = parser.parse_args()

    out_path = args.out or args.msa
    # write via a temporary file so a crash cannot truncate a finished alignment
    tmp_path = out_path + ".tmp"
    unmatched = fix(args.msa, args.pdb, tmp_path)
    os.replace(tmp_path, out_path)

    if unmatched:
        print(
            f"WARNING: {len(unmatched)} MSA header(s) match no PDB file in "
            f"{args.pdb}, left unchanged (e.g. {unmatched[:3]})",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
