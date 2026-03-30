import argparse
import os
import json
from typing import List, Set
import numpy as np
import pandas as pd
import glob


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
            print(tool)
            df = pd.read_csv(filename, index_col=False, sep=' ')
            df["tool"] = tool
            df["run_name"] = run_name
            #df = df.set_index(df.dataset + "_" + df.family)
            dfs.append(df)

    # Concatenate all dataframes into a single dataframe
    merged_df = pd.concat(dfs, ignore_index=True)

    return merged_df



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

    common_tools = find_tools(run_names, not args.all)

    if args.compare_tools:
        tool1, tool2 = args.compare_tools
        for t in [tool1, tool2]:
            if t not in common_tools:
                raise ValueError(f"Tool '{t}' not found in the given run(s). Available tools: {', '.join(sorted(common_tools))}")
        common_tools = {tool1, tool2}

    df = make_merged_df(run_names, common_tools)

    if args.compare_tools:
        tool1, tool2 = args.compare_tools

        df1 = df[df["tool"] == tool1].copy()
        df2 = df[df["tool"] == tool2].copy()

        merged = pd.merge(
            df1, df2,
            on=["run_name", "sample"],
            suffixes=(f"_{tool1}", f"_{tool2}")
        )

        merged["SP-Score_diff"] = merged[f"SP-Score_{tool2}"] - merged[f"SP-Score_{tool1}"]
        merged["TC_diff"] = merged[f"TC_{tool2}"] - merged[f"TC_{tool1}"]
        merged["runtime_diff"] = merged[f"s_{tool2}"] - merged[f"s_{tool1}"]

        merged_sorted = merged.sort_values(by="TC_diff", ascending=False)

        print(f"\nComparison: {tool2} vs {tool1}")
        print(f"Positive differences indicate {tool2} performed better (except runtime)\n")
        print(merged_sorted[["run_name", "sample", "SP-Score_diff", "TC_diff", "runtime_diff"]].to_string(index=False))

        print("\nSummary (mean differences by run):")
        print(merged.groupby(["run_name"])[["SP-Score_diff", "TC_diff", "runtime_diff"]].mean())

        if len(run_names) > 1:
            print("\nTotal (mean differences):")
            print(merged[["SP-Score_diff", "TC_diff", "runtime_diff"]].mean())

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
        merged["runtime_diff"] = merged[f"s_{run2_name}"] - merged[f"s_{run1_name}"]

        # Sort by TC difference (descending)
        merged_sorted = merged.sort_values(by="TC_diff", ascending=False)

        # Display results
        print(f"\nComparison: {run2_name} vs {run1_name}")
        print(f"Positive differences indicate {run2_name} performed better (except runtime)\n")
        print(merged_sorted[["tool", "sample", "SP-Score_diff", "TC_diff", "runtime_diff"]].to_string(index=False))

        # Print summary statistics
        print("\nSummary (mean differences by tool):")
        print(merged.groupby(["tool"])[["SP-Score_diff", "TC_diff", "runtime_diff"]].mean())

    elif args.detailed:
        df_sorted = df.sort_values(by="TC", ascending=False)
        print(df_sorted[["run_name", "tool", "sample", "SP-Score", "TC", "s", "h:m:s", "success"]].to_string(index=False))
    else:
        print(
            df.groupby(
                ["run_name", "tool"]
            )[["SP-Score", "TC", "s", "success"]].mean()
        )
        print("Total:")
        print(
            df.groupby(["tool"])[["SP-Score", "TC", "s", "success"]].mean()
        )
