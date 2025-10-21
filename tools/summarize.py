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


def find_tool_insersection(run_names: List[str]) -> Set[str]:
    """
    Find the common tools in the given run names.
    """
    tools: Set[str] | None = None
    for run in run_names:
        config_path = os.path.join("configs", run+".json")
        with open(config_path) as json_data:
            config = json.load(json_data)
            run_tools = config["tools"].keys()
            tools = tools.intersection(run_tools) if tools else set(run_tools)
    return tools if tools is not None else set()


def make_merged_df(run_names: List[str], tools: Set[str]) -> pd.DataFrame:
    # Create an empty list to store the dataframes
    dfs = []

    # Read each file and append its contents to the list
    for run_name in run_names:
        for i,tool in enumerate(tools):
            print(tool)
            filename = f"results/{run_name}/{tool}.out"
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
    
    print(f"Processing {len(run_names)} run(s): {', '.join(run_names)}")
    
    validate_run_names(run_names)

    common_tools = find_tool_insersection(run_names)

    df = make_merged_df(run_names, common_tools)

    if args.detailed:
        df_sorted = df.sort_values(by="TC", ascending=False)
        print(df_sorted[["run_name", "tool", "sample", "SP-Score", "TC", "s", "success"]].to_string(index=False))
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