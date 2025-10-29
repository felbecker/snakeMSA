import argparse
import glob
import os
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from learnMSA.msa_hmm.SequenceDataset import SequenceDataset

argparser = argparse.ArgumentParser(
    description="Show statistics for FASTA files in a directory."
)
argparser.add_argument(
    '-i', '--input', type=str, required=True,
    help="Directory containing FASTA files (supports wildcards with *)."
)
argparser.add_argument(
    '--pattern', type=str, default="*.fasta",
    help="Pattern to match FASTA files (default: *.fasta)."
)
argparser.add_argument(
    '--detailed', action='store_true',
    help="Show per-sequence statistics instead of just summaries."
)

args = argparser.parse_args()


def find_fasta_files(directory: str, pattern: str) -> List[Path]:
    """
    Find all FASTA files in the given directory matching the pattern.

    Args:
        directory: Directory to search for FASTA files.
        pattern: Glob pattern to match files (e.g., "*.fasta", "*.fa").

    Returns:
        List of Path objects for matching files.
    """
    search_path = os.path.join(directory, pattern)
    files = glob.glob(search_path)

    # Also try common FASTA extensions
    if not files:
        for ext in ["*.fa", "*.faa", "*.fna", "*.fasta"]:
            search_path = os.path.join(directory, ext)
            files.extend(glob.glob(search_path))

    return [Path(f) for f in sorted(set(files))]


def analyze_dataset(filepath: Path) -> dict:
    """
    Analyze a FASTA file and return statistics.

    Args:
        filepath: Path to the FASTA file.

    Returns:
        Dictionary containing statistics about the dataset.
    """
    try:
        dataset = SequenceDataset(filepath=filepath, fmt="fasta")

        stats = {
            "file": filepath.name,
            "num_sequences": len(dataset),
            "min_length": int(np.min(dataset.seq_lens)) if len(dataset) > 0 else 0,
            "avg_length": float(np.mean(dataset.seq_lens)) if len(dataset) > 0 else 0.0,
            "max_length": int(np.max(dataset.seq_lens)) if len(dataset) > 0 else 0,
            "total_residues": int(np.sum(dataset.seq_lens)),
            "std_length": float(np.std(dataset.seq_lens)) if len(dataset) > 0 else 0.0,
        }

        dataset.close()
        return stats

    except Exception as e:
        return {
            "file": filepath.name,
            "num_sequences": 0,
            "min_length": 0,
            "avg_length": 0.0,
            "max_length": 0,
            "total_residues": 0,
            "std_length": 0.0,
            "error": str(e),
        }


def make_statistics_df(files: List[Path]) -> pd.DataFrame:
    """
    Create a DataFrame with statistics for all FASTA files.

    Args:
        files: List of FASTA file paths.

    Returns:
        DataFrame with statistics for each file.
    """
    stats_list = []

    for filepath in files:
        print(f"Analyzing {filepath.name}...", file=sys.stderr)
        stats = analyze_dataset(filepath)
        stats_list.append(stats)

    df = pd.DataFrame(stats_list)
    return df


if __name__ == '__main__':
    # Find FASTA files
    if not os.path.exists(args.input):
        raise ValueError(f"Directory {args.input} does not exist.")

    fasta_files = find_fasta_files(args.input, args.pattern)

    if not fasta_files:
        raise ValueError(
            f"No FASTA files found in {args.input} matching pattern {args.pattern}"
        )

    print(f"Found {len(fasta_files)} FASTA file(s) in {args.input}", file=sys.stderr)
    print("", file=sys.stderr)

    # Analyze files
    df = make_statistics_df(fasta_files)

    # Display results
    if args.detailed:
        # Show all statistics for each file
        print(df.to_string(index=False))
    else:
        # Show summary statistics
        print(df[["file", "num_sequences", "min_length", "avg_length", "max_length"]].to_string(index=False))

        # Show overall statistics if multiple files
        if len(fasta_files) > 1:
            print("\n" + "="*80)
            print("OVERALL STATISTICS:")
            print("="*80)
            print(f"Total files:           {len(fasta_files)}")
            print(f"Total sequences:       {df['num_sequences'].sum()}")
            print(f"Total residues:        {df['total_residues'].sum()}")
            print(f"Average sequences/file: {df['num_sequences'].mean():.2f}")
            print(f"Min sequence length:   {df['min_length'].min()}")
            print(f"Max sequence length:   {df['max_length'].max()}")
            print(f"Overall avg length:    {df['avg_length'].mean():.2f}")
