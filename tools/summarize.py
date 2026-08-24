import argparse
import io
import os
import json
import re
import sys
from pathlib import Path
from typing import List, Set
import numpy as np
import pandas as pd
import glob

sys.path.insert(0, str(Path(__file__).parent))
from data_stats import compute_stats


STRUCT_SCORE_COLS = ["lddt", "lddt_core", "lddt_muscle", "dali_z"]

#: Structural scores reported unless --all-struct-scores is given.
DEFAULT_STRUCT_SCORE_COLS = ["lddt_muscle"]

#: Peak GPU memory in MiB, recorded by the align rule. -1 means "not
#: monitored" (CPU-only tool, or a run made before GPU tracking existed).
GPU_MEM_COL = "max_gpu_mem"


argparser = argparse.ArgumentParser(
    description="Summarize the results of one or more runs."
)
argparser.add_argument(
    '-i', nargs='+', type=str, required=True,
    help="One or more run names (supports wildcards with *)."
)
argparser.add_argument(
    '--detailed', action='store_true',
    help="Show all scores instead of just means."
)
argparser.add_argument(
    '--more_detailed', action='store_true',
    help="Show memory used as well as dataset properties."
)
argparser.add_argument(
    '--compare', action='store_true',
    help="Compare exactly two runs and show SP-Score and TC differences per sample."
)
argparser.add_argument(
    '--compare-tools', nargs=2, type=str, metavar=('TOOL1', 'TOOL2'),
    help="Compare exactly two tools within the same run(s) and show SP-Score "
    "and TC differences per sample."
)
argparser.add_argument(
    '--all', action='store_true',
    help="When multiple run names are given, shows all tools of all runs."\
    " Per default, only the intersection of tools is shown."
)
argparser.add_argument(
    '--tools', nargs='+', type=str, metavar='TOOL',
    help="Restrict output to this subset of tools (applies to all modes)."
)
argparser.add_argument(
    '--all-struct-scores', action='store_true', dest='all_struct_scores',
    help="Report all structural scores (lddt, lddt_core, lddt_muscle, dali_z) "
         "instead of only lddt_muscle."
)
argparser.add_argument(
    '--barplots', action='store_true',
    help="Create a bar plot (SP-Score, TC, runtime) for each run."
)
argparser.add_argument(
    '--barplot_lddt', nargs='?', const='lddt', default=None, metavar='LDDT_TYPE',
    help="Create a bar plot for LDDT. Optionally specify the variant to plot: "
         "lddt (default), lddt_core, lddt_muscle, or dali_z."
)
argparser.add_argument(
    '--average-replicas', action='store_true', dest='average_replicas',
    help="Average results across tools that share the same base name but differ only "
         "by a '_replica_N' suffix (e.g. tool_replica_1, tool_replica_2 → tool)."
)
argparser.add_argument(
    '--barplot_title', type=str, default=None, metavar='TITLE',
    help="Title to use for bar plots (overrides the default run name title)."
)
argparser.add_argument(
    '--barplot_fontsize', type=int, default=None, metavar='SIZE',
    help="Base font size for bar plots."
)
argparser.add_argument(
    '--barplot_tool_names', type=str, default=None, metavar='FILE',
    help="Tab-separated file mapping tool names: one row per tool with columns "
         "'name_in_dataframe<TAB>name_in_plot'. Whitespace is allowed in name_in_plot."
)

args = argparser.parse_args()

def expand_run_patterns(patterns: List[str]) -> List[str]:
    """
    Expand wildcard patterns in run names to actual run names.
    For example, 'benchmark_*' will expand to all runs starting with 'benchmark_'.
    """
    expanded_runs = []
    for pattern in patterns:
        if '*' in pattern:
            # Find matching config files
            config_pattern = os.path.join("configs", pattern + ".json")
            matching_configs = glob.glob(config_pattern)

            # Extract run names from config file paths
            for config_path in matching_configs:
                run_name = os.path.basename(config_path)[:-5]  # Remove .json
                # Check if corresponding results directory exists
                if os.path.exists(os.path.join("results", run_name)):
                    expanded_runs.append(run_name)
        else:
            expanded_runs.append(pattern)

    return expanded_runs

def validate_run_names(run_names: List[str]) -> None:
    """
    Validate the run names by checking if they exist.
    """
    for run in run_names:
        if (not os.path.exists("results/"+run)
            or not os.path.exists(os.path.join("configs", run+".json"))):
            raise ValueError(f"Run {run} does not exist.")


def find_tools(run_names: List[str], intersection: bool) -> Set[str]:
    """
    Find the common tools in the given run names.
    """
    tools: Set[str] | None = None
    for run in run_names:
        config_path = os.path.join("configs", run+".json")
        with open(config_path) as json_data:
            config = json.load(json_data)
            run_tools = config["tools"].keys()
            if intersection:
                tools = tools.intersection(run_tools) if tools else set(run_tools)
            else:
                tools = tools.union(run_tools) if tools else set(run_tools)
    return tools if tools is not None else set()


def make_merged_df(run_names: List[str], tools: Set[str]) -> pd.DataFrame:
    # Create an empty list to store the dataframes
    dfs = []

    # Read each file and append its contents to the list
    for run_name in run_names:
        for i,tool in enumerate(tools):
            filename = f"results/{run_name}/{tool}.out"
            if not os.path.isfile(filename):
                continue
            with open(filename) as fh:
                content = fh.read()
            # Normalize "N day[s], h:m:s" → "total_h:m:s" so the
            # whitespace-separated columns stay aligned, and replace tabs
            # with spaces (Snakemake benchmark files use tab separators).
            content = re.sub(
                r'(\d+) days?,\s*(\d+):(\d+):(\d+)',
                lambda m: f"{int(m.group(1)) * 24 + int(m.group(2))}:{m.group(3)}:{m.group(4)}",
                content
            )
            content = content.replace('\t', ' ')
            df = pd.read_csv(io.StringIO(content), index_col=False, sep=' ')
            df["tool"] = tool
            df["run_name"] = run_name
            #df = df.set_index(df.dataset + "_" + df.family)
            dfs.append(df)

    # Concatenate all dataframes into a single dataframe
    merged_df = pd.concat(dfs, ignore_index=True)

    # Treat -1 as "not available" so it is excluded from averages
    for _col in STRUCT_SCORE_COLS + [GPU_MEM_COL]:
        if _col in merged_df.columns:
            merged_df[_col] = merged_df[_col].replace(-1, np.nan)

    return merged_df



def wanted_struct_score_cols() -> List[str]:
    """The structural scores the user asked for: only lddt_muscle by default."""
    return STRUCT_SCORE_COLS if args.all_struct_scores else DEFAULT_STRUCT_SCORE_COLS


def struct_score_cols(frame: pd.DataFrame) -> List[str]:
    """The requested structural scores that are present in ``frame``."""
    return [c for c in wanted_struct_score_cols() if c in frame.columns]


def struct_diff_cols(frame: pd.DataFrame, suffix1: str, suffix2: str) -> List[str]:
    """The requested structural scores available on both sides of a comparison."""
    return [
        c for c in wanted_struct_score_cols()
        if f"{c}{suffix1}" in frame.columns and f"{c}{suffix2}" in frame.columns
    ]


def gpu_mem_cols(frame: pd.DataFrame) -> List[str]:
    """The peak GPU memory column, if any of the runs actually recorded it.

    Kept out of the output entirely for CPU-only runs and for older results
    files, where the column is either missing or all -1.
    """
    if GPU_MEM_COL in frame.columns and frame[GPU_MEM_COL].notna().any():
        return [GPU_MEM_COL]
    return []


if __name__ == '__main__':
    # Expand wildcard patterns in run names
    run_names = expand_run_patterns(args.i)

    if not run_names:
        raise ValueError("No matching runs found for the given patterns.")

    # Validate compare mode requires exactly 2 runs
    if args.compare and len(run_names) != 2:
        raise ValueError("--compare requires exactly 2 run names. Got {}".format(len(run_names)))

    if args.compare and args.compare_tools:
        raise ValueError("--compare and --compare-tools are mutually exclusive.")

    print(f"Processing {len(run_names)} run(s): {', '.join(run_names)}")

    validate_run_names(run_names)

    # Cross-run compare: exactly 2 runs + 2 tools → compare tool1 from run1 vs tool2 from run2
    cross_run_compare = bool(args.compare_tools) and len(run_names) == 2

    if cross_run_compare:
        tool1, tool2 = args.compare_tools
        run1_tools = find_tools([run_names[0]], intersection=True)
        run2_tools = find_tools([run_names[1]], intersection=True)
        if tool1 not in run1_tools:
            raise ValueError(
                f"Tool '{tool1}' not found in run '{run_names[0]}'. "
                f"Available: {', '.join(sorted(run1_tools))}"
            )
        if tool2 not in run2_tools:
            raise ValueError(
                f"Tool '{tool2}' not found in run '{run_names[1]}'. "
                f"Available: {', '.join(sorted(run2_tools))}"
            )
        common_tools = {tool1, tool2}
    else:
        common_tools = find_tools(run_names, not args.all)

        # Apply --tools filter if specified
        if args.tools:
            unknown = set(args.tools) - common_tools
            if unknown:
                print(f"Warning: the following tools were not found in the run(s) and will be ignored: {', '.join(sorted(unknown))}")
            common_tools = common_tools.intersection(set(args.tools))
            if not common_tools:
                raise ValueError("No matching tools remain after applying --tools filter.")

        if args.compare_tools:
            tool1, tool2 = args.compare_tools
            for t in [tool1, tool2]:
                if t not in common_tools:
                    raise ValueError(f"Tool '{t}' not found in the given run(s). Available tools: {', '.join(sorted(common_tools))}")
            common_tools = {tool1, tool2}

    df = make_merged_df(run_names, common_tools)

    if args.average_replicas:
        df["tool"] = df["tool"].str.replace(r'_replica_\d+$', '', regex=True)

    if args.compare_tools:
        tool1, tool2 = args.compare_tools

        if cross_run_compare:
            run1_name, run2_name = run_names[0], run_names[1]
            df1 = df[(df["tool"] == tool1) & (df["run_name"] == run1_name)].copy()
            df2 = df[(df["tool"] == tool2) & (df["run_name"] == run2_name)].copy()

            merged = pd.merge(df1, df2, on=["sample"], suffixes=("_1", "_2"))
            merged["SP-Score_diff"] = merged["SP-Score_2"] - merged["SP-Score_1"]
            merged["TC_diff"] = merged["TC_2"] - merged["TC_1"]
            _lddt_diff_cols = []
            for _lc in struct_diff_cols(merged, "_1", "_2"):
                merged[f"{_lc}_diff"] = merged[f"{_lc}_2"] - merged[f"{_lc}_1"]
                _lddt_diff_cols.append(f"{_lc}_diff")
            merged["runtime_diff"] = merged["s_2"] - merged["s_1"]

            label1 = f"{tool1} ({run1_name})"
            label2 = f"{tool2} ({run2_name})"
            print(f"\nComparison: {label2} vs {label1}")
            print(f"Positive differences indicate {label2} performed better (except runtime)\n")

            diff_cols = ["sample", "SP-Score_diff", "TC_diff"] + _lddt_diff_cols + ["runtime_diff"]

            if args.more_detailed:
                # Load dataset statistics from run1's config
                config_path = os.path.join("configs", run1_name + ".json")
                with open(config_path) as fh:
                    config = json.load(fh)
                data_path = Path(config.get("data_path", ""))
                tsv_unaligned = data_path / "statistics_unaligned.tsv"
                tsv_aligned = data_path / "statistics_aligned.tsv"
                if tsv_unaligned.exists():
                    stats_u = pd.read_csv(tsv_unaligned, sep="\t")
                else:
                    print(f"Computing unaligned statistics for {run1_name} ...", file=sys.stderr)
                    stats_u = compute_stats([str(data_path / "unaligned")], aligned=False)
                    stats_u.to_csv(tsv_unaligned, sep="\t", index=False)
                if tsv_aligned.exists():
                    stats_a = pd.read_csv(tsv_aligned, sep="\t")
                else:
                    print(f"Computing aligned statistics for {run1_name} ...", file=sys.stderr)
                    stats_a = compute_stats([str(data_path / "aligned")], aligned=True)
                    stats_a.to_csv(tsv_aligned, sep="\t", index=False)
                ds_stats = stats_u[["family", "num_seq", "max_len", "avg_len"]].merge(
                    stats_a[["family", "avg_pairwise_identity"]].rename(
                        columns={"avg_pairwise_identity": "pid"}
                    ),
                    on="family", how="outer",
                ).rename(columns={"family": "sample"})
                merged = merged.merge(ds_stats, on="sample", how="left")
                diff_cols = ["sample", "num_seq", "max_len", "avg_len", "pid",
                             "SP-Score_diff", "TC_diff"] + _lddt_diff_cols + ["runtime_diff"]

            if args.detailed or args.more_detailed:
                merged_sorted = merged.sort_values(by="TC_diff", ascending=False)
                print(merged_sorted[diff_cols].to_string(index=False))

            print("\nSummary (mean differences):")
            _summary_diff_cols = ["SP-Score_diff", "TC_diff"] + _lddt_diff_cols + ["runtime_diff"]
            print(merged[_summary_diff_cols].mean().to_string())

        else:
            df1 = df[df["tool"] == tool1].copy()
            df2 = df[df["tool"] == tool2].copy()

            merged = pd.merge(
                df1, df2,
                on=["run_name", "sample"],
                suffixes=(f"_{tool1}", f"_{tool2}")
            )

            merged["SP-Score_diff"] = merged[f"SP-Score_{tool2}"] - merged[f"SP-Score_{tool1}"]
            merged["TC_diff"] = merged[f"TC_{tool2}"] - merged[f"TC_{tool1}"]
            _lddt_diff_cols2 = []
            for _lc in struct_diff_cols(merged, f"_{tool1}", f"_{tool2}"):
                merged[f"{_lc}_diff"] = merged[f"{_lc}_{tool2}"] - merged[f"{_lc}_{tool1}"]
                _lddt_diff_cols2.append(f"{_lc}_diff")
            merged["runtime_diff"] = merged[f"s_{tool2}"] - merged[f"s_{tool1}"]

            merged_sorted = merged.sort_values(by="TC_diff", ascending=False)

            print(f"\nComparison: {tool2} vs {tool1}")
            print(f"Positive differences indicate {tool2} performed better (except runtime)\n")
            print(merged_sorted[["run_name", "sample", "SP-Score_diff", "TC_diff"] + _lddt_diff_cols2 + ["runtime_diff"]].to_string(index=False))

            print("\nSummary (mean differences by run):")
            print(merged.groupby(["run_name"])[["SP-Score_diff", "TC_diff"] + _lddt_diff_cols2 + ["runtime_diff"]].mean())

            if len(run_names) > 1:
                print("\nTotal (mean differences):")
                print(merged[["SP-Score_diff", "TC_diff"] + _lddt_diff_cols2 + ["runtime_diff"]].mean())

    elif args.compare:
        # Compare mode: show differences between two runs
        run1_name, run2_name = run_names[0], run_names[1]

        # Separate data for each run
        df1 = df[df["run_name"] == run1_name].copy()
        df2 = df[df["run_name"] == run2_name].copy()

        # Merge on tool and sample to compute differences
        merged = pd.merge(
            df1, df2,
            on=["tool", "sample"],
            suffixes=(f"_{run1_name}", f"_{run2_name}")
        )

        # Compute differences (run2 - run1)
        merged["SP-Score_diff"] = merged[f"SP-Score_{run2_name}"] - merged[f"SP-Score_{run1_name}"]
        merged["TC_diff"] = merged[f"TC_{run2_name}"] - merged[f"TC_{run1_name}"]
        _lddt_diff_cols3 = []
        for _lc in struct_diff_cols(merged, f"_{run1_name}", f"_{run2_name}"):
            merged[f"{_lc}_diff"] = merged[f"{_lc}_{run2_name}"] - merged[f"{_lc}_{run1_name}"]
            _lddt_diff_cols3.append(f"{_lc}_diff")
        merged["runtime_diff"] = merged[f"s_{run2_name}"] - merged[f"s_{run1_name}"]

        # Sort by TC difference (descending)
        merged_sorted = merged.sort_values(by="TC_diff", ascending=False)

        # Display results
        print(f"\nComparison: {run2_name} vs {run1_name}")
        print(f"Positive differences indicate {run2_name} performed better (except runtime)\n")
        print(merged_sorted[["tool", "sample", "SP-Score_diff", "TC_diff"] + _lddt_diff_cols3 + ["runtime_diff"]].to_string(index=False))

        # Print summary statistics
        print("\nSummary (mean differences by tool):")
        print(merged.groupby(["tool"])[["SP-Score_diff", "TC_diff"] + _lddt_diff_cols3 + ["runtime_diff"]].mean())

    elif args.detailed or args.more_detailed:
        df_sorted = df.sort_values(by="TC", ascending=False)
        if args.more_detailed:
            # Load or compute per-family dataset statistics and merge into df
            stats_parts = []
            for run in run_names:
                config_path = os.path.join("configs", run + ".json")
                with open(config_path) as fh:
                    config = json.load(fh)
                data_path = Path(config.get("data_path", ""))
                tsv_unaligned = data_path / "statistics_unaligned.tsv"
                tsv_aligned = data_path / "statistics_aligned.tsv"
                if not tsv_unaligned.exists():
                    print(
                        f"Computing unaligned statistics for {run} ...",
                        file=sys.stderr,
                    )
                    unaligned_dir = str(data_path / "unaligned")
                    stats_u = compute_stats([unaligned_dir], aligned=False)
                    stats_u.to_csv(tsv_unaligned, sep="\t", index=False)
                else:
                    stats_u = pd.read_csv(tsv_unaligned, sep="\t")
                if not tsv_aligned.exists():
                    print(
                        f"Computing aligned statistics for {run} ...",
                        file=sys.stderr,
                    )
                    aligned_dir = str(data_path / "aligned")
                    stats_a = compute_stats([aligned_dir], aligned=True)
                    stats_a.to_csv(tsv_aligned, sep="\t", index=False)
                else:
                    stats_a = pd.read_csv(tsv_aligned, sep="\t")
                # Combine the columns we want from both TSVs
                ds_stats = stats_u[["family", "num_seq", "max_len", "avg_len"]].merge(
                    stats_a[["family", "avg_pairwise_identity"]],
                    on="family",
                    how="outer",
                )
                ds_stats = ds_stats.rename(columns={"family": "sample"})
                ds_stats["run_name"] = run
                stats_parts.append(ds_stats)
            all_stats = pd.concat(stats_parts, ignore_index=True)
            df_sorted = df_sorted.merge(
                all_stats, on=["run_name", "sample"], how="left"
            )
            df_sorted.rename(
                columns={"avg_pairwise_identity": "pid"}, inplace=True
            )
            _lddt_cols = struct_score_cols(df_sorted)
            print(df_sorted[
                ["run_name", "tool", "sample",
                 "num_seq", "max_len", "avg_len", "pid",
                 "SP-Score", "TC"] + _lddt_cols + ["s", "h:m:s", "success"]
            ].to_string(index=False))
        else:
            _lddt_cols = struct_score_cols(df_sorted)
            _gpu_cols = gpu_mem_cols(df_sorted)
            print(df_sorted[
                ["run_name", "tool", "sample", "SP-Score", "TC"] + _lddt_cols
                + ["s", "h:m:s", "max_rss"] + _gpu_cols + ["success"]
            ].to_string(index=False))
    else:
        with pd.option_context('display.max_columns', None, 'display.width', None):
            _lddt_cols = struct_score_cols(df)
            _cols = ["SP-Score", "TC"] + _lddt_cols + ["s"] + gpu_mem_cols(df) + ["success"]
            print(
                df.groupby(
                    ["run_name", "tool"]
                )[_cols].mean()
            )
            print("Total:")
            print(
                df.groupby(["tool"])[_cols].mean()
            )

    if args.barplots or args.barplot_lddt:
        from plots import barplot, barplot_lddt
        tool_order = args.tools if args.tools else None
        if len(run_names) == 1:
            name = run_names[0]
        else:
            name = "_".join(run_names)
        # Parse optional tool name mapping file
        tool_name_map = None
        if args.barplot_tool_names:
            tool_name_map = {}
            with open(args.barplot_tool_names) as _fh:
                for _line in _fh:
                    _line = _line.rstrip("\n")
                    if not _line or _line.startswith("#"):
                        continue
                    _parts = _line.split(None, 1)
                    if len(_parts) == 2:
                        tool_name_map[_parts[0]] = _parts[1]
        plot_title = args.barplot_title if args.barplot_title else name
        plot_fontsize = args.barplot_fontsize
        if args.barplots:
            barplot(df, run_name=name, tools=tool_order, output_path=f"plots/{name}_barplot.pdf",
                    title=plot_title, font_size=plot_fontsize, tool_name_map=tool_name_map)
        if args.barplot_lddt:
            lddt_col = args.barplot_lddt
            barplot_lddt(df, run_name=name, tools=tool_order, lddt_col=lddt_col,
                         output_path=f"plots/{name}_barplot_{lddt_col}.pdf",
                         title=plot_title, font_size=plot_fontsize, tool_name_map=tool_name_map)
