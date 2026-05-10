import argparse
import sys
from pathlib import Path

import pandas as pd

from learnMSA.util import SequenceDataset, AlignedDataset


def collect_families(inputs: list[str]) -> dict[str, str]:
    """Expand a list of files/directories into a {family_name: path} dict."""
    families: dict[str, str] = {}
    for p in inputs:
        path = Path(p)
        if path.is_dir():
            for f in sorted(path.iterdir()):
                if f.is_file():
                    families[f.stem] = str(f)
        else:
            families[path.stem] = str(path)
    return families


def compute_stats(inputs: list[str], aligned: bool) -> pd.DataFrame:
    """
    Compute sequence statistics for the given input files/directories.

    Args:
        inputs: List of file paths or directories containing FASTA files.
        aligned: If True, treat as aligned MSAs and also compute average
            pairwise identity.

    Returns:
        DataFrame with columns: family, num_seq, min_len, max_len, avg_len
        and, when aligned=True, avg_pairwise_identity.
    """
    families = collect_families(inputs)
    rows = []
    for family, path in families.items():
        if aligned:
            ds = AlignedDataset(filepath=path, validate_alphabet=False)
        else:
            ds = SequenceDataset(filepath=path, validate_alphabet=False)
        ds.validate_dataset(single_seq_ok=True)
        row: dict = {
            "family": family,
            "num_seq": ds.num_seq,
            "min_len": int(ds.seq_lens.min()),
            "max_len": int(ds.seq_lens.max()),
            "avg_len": round(float(ds.seq_lens.mean()), 2),
        }
        if aligned:
            row["avg_pairwise_identity"] = round(
                ds.average_pairwise_identity(), 4
            )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report sequence statistics per family (fasta file)."
    )
    parser.add_argument(
        "-i", "--input",
        nargs="+",
        metavar="FILE",
        required=True,
        help="One or more fasta files or directories. The file basename "
             "(without extension) is used as the family name.",
    )
    parser.add_argument(
        "--aligned",
        action="store_true",
        help="Treat inputs as multiple sequence alignments and also report "
             "average pairwise identity.",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="FILE",
        default=None,
        help="Write TSV output to FILE instead of stdout.",
    )
    args = parser.parse_args()

    df = compute_stats(args.input, args.aligned)

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        df.to_csv(out, sep="\t", index=False)
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()


