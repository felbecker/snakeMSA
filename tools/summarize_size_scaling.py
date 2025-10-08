import argparse
import os
import re
from typing import Tuple, Optional
import numpy as np
import pandas as pd


argparser = argparse.ArgumentParser(
    description="Summarize the results of a size scaling experiment."
)
argparser.add_argument(
    '-i', nargs='+', type=str, required=True,
    help="Path(s) to one or more input data files."
)
argparser.add_argument(
    '-o', '--output', type=str, required=False,
    help="Optional output file path (only used for single input)."
)
argparser.add_argument(
    '--detailed', action='store_true',
    help="Show detailed breakdown by family and sequence count."
)

args = argparser.parse_args()


def parse_sample_key(sample: str) -> Tuple[str, int, int]:
    """
    Parse the sample key to extract family name, sequence count, and run number.
    Example: "lyase_1_500_1" -> family="lyase_1", seq_count=500, run=1
    """
    # Split by underscore and find the pattern where we have digits followed by more digits
    parts = sample.split('_')

    # The last part should be the run number
    run_num = int(parts[-1])

    # The second to last should be the sequence count
    seq_count = int(parts[-2])

    # Everything before that is the family name
    family = '_'.join(parts[:-2])

    return family, seq_count, run_num


def load_and_process_data(input_file: str) -> pd.DataFrame:
    """
    Load the data file and add parsed columns for family, seq_count, and run_num.
    """
    # Read the data file
    df = pd.read_csv(input_file, sep=' ')

    # Parse the sample column
    parsed = df['sample'].apply(parse_sample_key)
    df['family'] = parsed.apply(lambda x: x[0])
    df['seq_count'] = parsed.apply(lambda x: x[1])
    df['run_num'] = parsed.apply(lambda x: x[2])

    return df


def summarize_by_seq_count(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute means over families and runs for each sequence count.
    """
    # Columns to display in output
    display_cols = ['SP-Score', 'Modeler', 'TC', 's']

    # Filter to only existing columns
    cols_to_aggregate = [col for col in display_cols if col in df.columns]

    # Group by sequence count and compute means
    summary = df.groupby('seq_count')[cols_to_aggregate].mean()

    # Sort by sequence count
    summary = summary.sort_index()

    return summary


def print_summary(
        df: pd.DataFrame,
        filename: Optional[str] = None,
        show_detailed: bool = False
) -> None:
    """
    Print summary statistics.

    Args:
        df: DataFrame with the data
        filename: Optional filename to display in header
        show_detailed: If True, show detailed breakdown by family
    """
    if filename:
        print("\n" + "="*80)
        print(f"FILE: {filename}")

    print("="*80 + "\n")

    summary = summarize_by_seq_count(df)
    print(summary)

    if show_detailed:
        print("\n" + "="*80)
        print("DETAILED BREAKDOWN BY FAMILY")
        print("="*80 + "\n")

        # Show breakdown by family and sequence count
        display_cols = ['SP-Score', 'Modeler', 'TC', 's']
        cols_to_show = [col for col in display_cols if col in df.columns]

        family_summary = df.groupby(['family', 'seq_count'])[cols_to_show].mean()
        print(family_summary)


if __name__ == '__main__':
    # Process each input file
    for i, input_file in enumerate(args.i):
        # Load and process the data
        df = load_and_process_data(input_file)

        # Print summary (show filename if multiple files)
        filename = input_file if len(args.i) > 1 else None
        print_summary(df, filename=filename, show_detailed=args.detailed)

        # Save to output file if specified (only for single input)
        if args.output and len(args.i) == 1:
            summary = summarize_by_seq_count(df)
            summary.to_csv(args.output)
            print(f"\nSummary saved to: {args.output}")
        elif args.output and len(args.i) > 1:
            # For multiple files, append index to output filename
            base, ext = os.path.splitext(args.output)
            output_file = f"{base}_{i}{ext}"
            summary = summarize_by_seq_count(df)
            summary.to_csv(output_file)
            print(f"\nSummary saved to: {output_file}")

