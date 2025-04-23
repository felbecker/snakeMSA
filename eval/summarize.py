import argparse 
import os
import json
import numpy as np
import pandas as pd


argparser = argparse.ArgumentParser(description='Summarize the results of one or more runs.')
argparser.add_argument('-i', nargs='+', type=str, required=True, help='One or more run names.')

args = argparser.parse_args()



def validate_run_names(run_names):
    """
    Validate the run names by checking if they exist.
    """
    for run in run_names:
        if (not os.path.exists("results/"+run) 
            or not os.path.exists(os.path.join("configs", run+".json"))):
            raise ValueError(f"Run {run} does not exist.")


def find_tool_insersection(run_names):
    """
    Find the common tools in the given run names.
    """
    tools = None
    for run in run_names:
        config_path = os.path.join("configs", run+".json")
        with open(config_path) as json_data:
            config = json.load(json_data)
            run_tools = config["tools"].keys()
            tools = tools.intersection(run_tools) if tools else set(run_tools)
    return tools


def make_merged_df(run_names, tools):
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
    validate_run_names(args.i)

    common_tools = find_tool_insersection(args.i)

    df = make_merged_df(args.i, common_tools)

    print(df.groupby(["run_name", "tool"])[["SP-Score", "TC", "s"]].mean())