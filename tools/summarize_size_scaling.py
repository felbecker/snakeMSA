"""Summarize a size scaling benchmark per alignment depth.

Averages the scores of every depth (number of aligned sequences) over all
families, subsamples and replicas of each tool. The samples of the benchmark
are named <family>_<depth>_<subsample>, e.g. rvp_200000_3. Prints a table
per tool and optionally plots a score over the depth with one curve per tool.

    python tools/summarize_size_scaling.py --plot --score TC
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plots import depth_curve
from summarize_missing_struct import read_tool, read_tool_names

SCORE_COLS = ["SP-Score", "Modeler", "TC", "lddt", "lddt_muscle",
              "lddt_muscle_L", "lddt_muscle_conf", "lddt_muscle_L_conf",
              "coverage", "lddt_fm", "lddt_ref", "lddt_core", "dali_z"]
RUNTIME_COL = "tool_s"


def read_scaling_tool(run: str, tool: str) -> pd.DataFrame:
    """The results of tool in run with the parsed columns ``family``,
    ``depth`` and ``subsample`` and the tool name in ``tool``."""
    df = read_tool(run, tool)
    if df is None:
        raise SystemExit(f"No results of {tool} in results/{run}.")
    parts = df["sample"].str.extract(r"^(.+)_(\d+)_(\d+)$")
    if parts.isna().any(axis=None):
        bad = df["sample"][parts.isna().any(axis=1)].iloc[0]
        raise SystemExit(f"Sample {bad} in {run} is not named "
                         "<family>_<depth>_<subsample>.")
    df["family"] = parts[0]
    df["depth"] = parts[1].astype(int)
    df["subsample"] = parts[2].astype(int)
    df["tool"] = tool
    # -1 marks failed alignments or missing scores
    for col in SCORE_COLS[:3] + [RUNTIME_COL]:
        if col in df.columns:
            df[col] = df[col].mask(df[col] < 0)
    return df


def summary_table(df: pd.DataFrame, cols: list, score: str,
                  by: list) -> pd.DataFrame:
    """Means of cols grouped by the columns by, with the number of
    alignments, of alignments with a value of score and of replicas, and the
    standard deviation of the per-replica means of score."""
    groups = df.groupby(by)
    table = groups[cols].mean()
    table.insert(0, "n", groups.size())
    table.insert(1, f"n_{score}", groups[score].count())
    table.insert(2, "replicas", groups["replica"].nunique())
    replica_means = df.groupby(by + ["replica"])[score].mean()
    table[f"{score}_sd"] = replica_means.groupby(level=by).std(ddof=1)
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--run", default="learnMSA_scale",
        help="Benchmark run with samples <family>_<depth>_<subsample>.")
    parser.add_argument(
        "--tools", nargs="+", default=["learnMSA", "learnMSA_struct"],
        help="Tools without their '_replica_N' suffix; replicas are "
             "averaged.")
    parser.add_argument(
        "--score", default="TC", choices=SCORE_COLS + [RUNTIME_COL],
        help="Score on the y axis of the plot and of the standard deviation "
             "column of the table.")
    parser.add_argument(
        "--detailed", action="store_true",
        help="Also print the means per family and depth.")
    parser.add_argument(
        "--csv", type=str, default=None, metavar="PATH",
        help="Save the table of all tools as CSV.")
    parser.add_argument(
        "--plot", action="store_true",
        help="Plot --score over the alignment depth, one curve per tool.")
    parser.add_argument(
        "--tool_names", type=str, default=None, metavar="PATH",
        help="File with display names for the legend, as in summarize.py "
             "--barplot_tool_names (e.g. tools/tool_names.txt).")
    parser.add_argument(
        "--title", type=str, default=None, metavar="TITLE",
        help="Title of the plot (default: none).")
    parser.add_argument(
        "--fontsize", type=int, default=None, metavar="SIZE",
        help="Base font size of the plot.")
    parser.add_argument(
        "--logy", action=argparse.BooleanOptionalAction, default=None,
        help="Logarithmic y axis (default: only for the runtime).")
    parser.add_argument(
        "--output", type=str, default=None, metavar="PATH",
        help="Where to save the plot. A path with an extension (e.g. "
             "figs/my_plot.pdf) is used as the file name; anything else is "
             "treated as a directory that receives the default file name. "
             "Default: plots/.")
    args = parser.parse_args()

    df = pd.concat([read_scaling_tool(args.run, tool) for tool in args.tools],
                   ignore_index=True)
    cols = [col for col in SCORE_COLS + [RUNTIME_COL] if col in df.columns]
    if args.score not in df.columns:
        raise SystemExit(f"No column {args.score} in the results.")

    failed = df[(df["success"] != 1) | df[["SP-Score", "TC"]].isna().any(
        axis=1)]
    if len(failed):
        print(f"Warning: {len(failed)} alignments failed, e.g. "
              f"{failed['tool'].iloc[0]} {failed['sample'].iloc[0]}.",
              file=sys.stderr)

    table = summary_table(df, cols, args.score, ["tool", "depth"])
    for tool in args.tools:
        sub = table.loc[tool]
        if sub["n"].nunique() > 1:
            print(f"Warning: the number of alignments of {tool} differs "
                  "between the depths.", file=sys.stderr)
        print(f"\n{tool} by alignment depth "
              f"({df[df['tool'] == tool]['family'].nunique()} families):")
        print(sub.to_string(float_format=lambda x: f"{x:.4f}"))
        if args.detailed:
            detail = summary_table(df[df["tool"] == tool], cols, args.score,
                                   ["family", "depth"])
            print(f"\n{tool} by family and depth:")
            print(detail.to_string(float_format=lambda x: f"{x:.4f}"))

    if args.csv is not None:
        os.makedirs(os.path.dirname(args.csv) or ".", exist_ok=True)
        table.to_csv(args.csv)
        print(f"\nSaved the table to {args.csv}")

    if args.plot:
        name = f"{args.run}_{args.score}.pdf"
        out = args.output
        if out is None:
            output_path = os.path.join("plots", name)
        elif os.path.splitext(out)[1]:
            output_path = out
        else:
            output_path = os.path.join(out, name)
        logy = args.logy if args.logy is not None \
            else args.score == RUNTIME_COL
        depth_curve(df, score_col=args.score, output_path=output_path,
                    tools=args.tools, title=args.title,
                    font_size=args.fontsize,
                    tool_name_map=read_tool_names(args.tool_names),
                    logy=logy)


if __name__ == "__main__":
    main()
