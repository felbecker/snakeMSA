"""Summarize the missing structure benchmark per fraction p of structures.

Averages the scores of every p over all clusters and replicas. The samples
of the benchmark run are named <cluster>_p<100p> (see make_missing_struct.py);
the base run with all structures (p = 1) uses the cluster names. Prints a
table and plots a score over p, by default together with other tools on the
same clusters as single points: tools that use structures at p = 1, tools
that use only sequences at p = 0.

    python tools/summarize_missing_struct.py
"""

import argparse
import glob
import io
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plots import missing_struct_curve

SCORE_COLS = ["lddt", "lddt_core", "lddt_muscle", "dali_z"]


def read_out(path: str) -> pd.DataFrame:
    """A results file of snakeMSA as a dataframe."""
    with open(path) as f:
        content = f.read()
    # "N day[s], h:m:s" -> "total_h:m:s", as in summarize.py
    content = re.sub(
        r"(\d+) days?,\s*(\d+):(\d+):(\d+)",
        lambda m: f"{int(m[1]) * 24 + int(m[2])}:{m[3]}:{m[4]}",
        content,
    ).replace("\t", " ")
    df = pd.read_csv(io.StringIO(content), index_col=False, sep=" ")
    # -1 means "not available"
    df[SCORE_COLS] = df[SCORE_COLS].replace(-1, np.nan)
    return df


def read_tool(run: str, tool: str) -> pd.DataFrame | None:
    """The results of tool in run, one row per alignment, with the replica
    number in the column ``replica`` (0 for a tool without replicas). The
    results of all replicas tool_replica_N are combined. None if run has
    no results of tool."""
    frames = []
    path = f"results/{run}/{tool}.out"
    if os.path.isfile(path):
        frames.append(read_out(path).assign(replica=0))
    pattern = re.compile(rf"{re.escape(tool)}_replica_(\d+)\.out")
    for path in sorted(glob.glob(f"results/{run}/{tool}_replica_*.out")):
        match = pattern.fullmatch(os.path.basename(path))
        if match is not None:
            frames.append(read_out(path).assign(replica=int(match[1])))
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def uses_structures(run: str, tool: str) -> bool:
    """Whether tool gets structures in the snakeMSA config of run."""
    with open(f"configs/{run}.json") as f:
        tools = json.load(f)["tools"]
    params = tools.get(tool) or tools.get(f"{tool}_replica_1")
    if params is None:
        raise SystemExit(f"configs/{run}.json has no tool {tool}.")
    params = " ".join(params["params"])
    return "{pdb}" in params or "--struct" in params


def read_tool_names(path: str | None) -> dict:
    """Display names as in summarize.py --barplot_tool_names."""
    names = {}
    if path is None:
        return names
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                names[parts[0]] = parts[1].replace("\\n", "\n")
    return names


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--run", default="learnMSA_pub_afdb_clusters_missing_struct",
        help="Benchmark run with samples <cluster>_p<100p>.")
    parser.add_argument(
        "--base_run", default="learnMSA_pub_afdb_clusters",
        help="Run with all structures (p = 1).")
    parser.add_argument(
        "--tool", default="learnMSA_struct",
        help="Tool name without its '_replica_N' suffix.")
    parser.add_argument(
        "--score", default="lddt_muscle", choices=SCORE_COLS,
        help="Score on the y axis of the plot.")
    parser.add_argument(
        "--compare", action=argparse.BooleanOptionalAction, default=True,
        help="Also draw the tools of --compare_tools as single points, with "
             "a legend: at p = 1 if they use structures, else at p = 0.")
    parser.add_argument(
        "--compare_runs", nargs="+",
        default=["learnMSA_pub_afdb_clusters",
                 "competition_pub_afdb_clusters"],
        help="Runs with the results of the --compare_tools.")
    parser.add_argument(
        "--compare_tools", nargs="+",
        default=["learnMSA", "famsa", "foldmason", "muscle_3d",
                 "muscle_3d_super7"],
        help="Tools to compare with, without '_replica_N' suffix; replicas "
             "are averaged.")
    parser.add_argument(
        "--tool_names", type=str, default=None, metavar="PATH",
        help="File with display names for the legend, as in summarize.py "
             "--barplot_tool_names (e.g. tools/tool_names_pdb.txt).")
    parser.add_argument(
        "--label", type=str, default=None, metavar="LABEL",
        help="Legend label of the curve (default: the display name of "
             "--tool).")
    parser.add_argument(
        "--title", type=str, default=None, metavar="TITLE",
        help="Title of the plot (default: none).")
    parser.add_argument(
        "--fontsize", type=int, default=None, metavar="SIZE",
        help="Base font size of the plot.")
    parser.add_argument(
        "--output", type=str, default=None, metavar="PATH",
        help="Where to save the plot. A path with an extension (e.g. "
             "figs/my_plot.pdf) is used as the file name; anything else is "
             "treated as a directory that receives the default file name. "
             "Default: plots/.")
    args = parser.parse_args()

    df = read_tool(args.run, args.tool)
    base = read_tool(args.base_run, args.tool)
    for frame, run in ((df, args.run), (base, args.base_run)):
        if frame is None:
            raise SystemExit(f"No results of {args.tool} in results/{run}.")
    suffix = df["sample"].str.extract(r"_p(\d+)$")[0]
    if suffix.isna().any():
        raise SystemExit(
            f"Samples without a _p<100p> suffix in {args.run}, e.g. "
            f"{df['sample'][suffix.isna()].iloc[0]}.")
    df["p"] = suffix.astype(int) / 100
    base["p"] = 1.0
    df = pd.concat([base, df], ignore_index=True)

    others = []
    if args.compare:
        for tool in args.compare_tools:
            found = [(run, read_tool(run, tool)) for run in args.compare_runs]
            found = [(run, frame) for run, frame in found if frame is not None]
            if not found:
                raise SystemExit(f"No results of {tool} in the runs "
                                 f"{' '.join(args.compare_runs)}.")
            run, frame = found[0]
            frame["tool"] = tool
            frame["p"] = 1.0 if uses_structures(run, tool) else 0.0
            others.append(frame)
    others = pd.concat(others, ignore_index=True) if others else None

    for frame in (df, others):
        if frame is None:
            continue
        failed = frame[(frame["success"] != 1)
                       | frame[SCORE_COLS].isna().any(axis=1)]
        if len(failed):
            print(f"Warning: {len(failed)} alignments failed or lack scores, "
                  f"e.g. {failed['sample'].iloc[0]}.", file=sys.stderr)

    groups = df.groupby("p")
    table = groups[SCORE_COLS].mean()
    table.insert(0, "n", groups.size())
    table.insert(1, "replicas", groups["replica"].nunique())
    # Error of the replicas as in the plot: the standard deviation of the
    # per-replica means.
    replica_means = df.groupby(["p", "replica"])[args.score].mean()
    table[f"{args.score}_sd"] = replica_means.groupby(level="p").std(ddof=1)
    table = table.sort_index(ascending=False)
    if table["n"].nunique() > 1:
        print("Warning: the number of alignments differs between the "
              "values of p.", file=sys.stderr)
    print(f"{args.tool} by fraction p of available structures:")
    print(table.to_string(float_format=lambda x: f"{x:.4f}"))

    if others is not None:
        groups = others.groupby("tool", sort=False)
        table = groups[SCORE_COLS].mean()
        table.insert(0, "p", groups["p"].first())
        table.insert(1, "n", groups.size())
        table.insert(2, "replicas", groups["replica"].nunique())
        print("\nOther tools (drawn at p = 1 if they use structures):")
        print(table.to_string(float_format=lambda x: f"{x:.4f}"))

    name = f"{args.run}_{args.score}.pdf"
    out = args.output
    if out is None:
        output_path = os.path.join("plots", name)
    elif os.path.splitext(out)[1]:
        output_path = out
    else:
        output_path = os.path.join(out, name)
    missing_struct_curve(df, score_col=args.score, output_path=output_path,
                         title=args.title, font_size=args.fontsize,
                         tool=args.tool, others=others,
                         tool_name_map=read_tool_names(args.tool_names),
                         curve_label=args.label)


if __name__ == "__main__":
    main()
