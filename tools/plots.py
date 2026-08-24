"""Bar plots for benchmark results of a single run with multiple tools."""

import os
from typing import List, Optional

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def barplot(
    df: pd.DataFrame,
    run_name: str,
    tools: Optional[List[str]] = None,
    output_path: Optional[str] = None,
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool_name_map: Optional[dict] = None,
) -> None:
    """Create a bar plot for a single run overlaying SP-Score (faded) and TC
    (solid) bars, with mean runtime as a label over each bar group.

    Parameters
    ----------
    df:
        Merged dataframe as returned by ``make_merged_df`` from summarize.py,
        already filtered to this run.
    run_name:
        Name used for the default output filename.
    tools:
        Ordered list of tool names to include. If *None* all tools in *df* are
        used, sorted by mean TC score descending.
    output_path:
        File path to save the figure (PDF). If *None* the figure is saved as
        ``<run_name>_barplot.pdf`` in the current working directory.
    title:
        Title for the figure. Defaults to *run_name*.
    font_size:
        Base font size. Defaults to 10.
    tool_name_map:
        Mapping from dataframe tool names to display names shown on the plot.
    """
    sns.set_style("whitegrid")

    available = df["tool"].unique().tolist()

    if tools is None:
        candidates = available
    else:
        candidates = [t for t in tools if t in available]

    # Always sort by mean TC descending
    order = (
        df[df["tool"].isin(candidates)].groupby("tool")["TC"].mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    if not order:
        raise ValueError("No matching tools found in the dataframe.")

    data = df[df["tool"].isin(order)].copy()
    means = data.groupby("tool")[["SP-Score", "TC", "s"]].mean()

    # Apply display-name mapping
    if tool_name_map:
        display_order = [tool_name_map.get(t, t) for t in order]
        data["tool"] = data["tool"].map(lambda t: tool_name_map.get(t, t))
        means.index = [tool_name_map.get(t, t) for t in means.index]
    else:
        display_order = order

    n_tools = len(display_order)
    font_size = font_size if font_size is not None else 10
    fig, ax = plt.subplots(figsize=(max(5, n_tools * 1.2), 4.5))
    fig.suptitle(title if title is not None else run_name, fontsize=font_size + 2, fontweight="bold")

    # Assign colors by stable alphabetical order so the same tool always
    # gets the same color regardless of the display ordering.
    stable_order = sorted(display_order)
    palette_deep = sns.color_palette("deep", n_tools)
    palette_dark = sns.color_palette("dark", n_tools)
    color_deep = {t: palette_deep[i] for i, t in enumerate(stable_order)}
    color_dark = {t: palette_dark[i] for i, t in enumerate(stable_order)}

    # SP-Score bars (faded)
    sns.barplot(
        data=data, x="tool", y="SP-Score", hue="tool", order=display_order,
        err_kws={"color": ".4", "linewidth": 1.5}, capsize=0.2,
        linewidth=1, edgecolor=".6",
        palette=color_deep, alpha=0.5, ax=ax,
    )

    # TC bars (solid) overlaid
    sns.barplot(
        data=data, x="tool", y="TC", hue="tool", order=display_order,
        errorbar=None,
        palette=color_dark, alpha=0.5, ax=ax,
    )

    # Runtime labels centred on each bar group
    y_min = means["TC"].min()
    y_mid = means["SP-Score"].max()
    label_y = y_min + (y_mid - y_min) / 2

    for tool in display_order:
        seconds = means.loc[tool, "s"]
        hours = seconds / 3600
        if hours >= 1:
            label = f"{hours:.2f} h"
        elif seconds >= 60:
            label = f"{seconds / 60:.1f} m"
        else:
            label = f"{seconds:.0f} s"
        ax.text(
            tool, label_y, label,
            ha="center", va="center",
            color="black", fontsize=font_size * 0.85,
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
        )

    ax.set_xlabel("")
    ax.set_ylabel("SP-Score (fade) / TC (solid)", fontsize=font_size)
    ax.set_ylim(bottom=max(0, y_min - 0.05))
    for item in ax.get_xticklabels():
        item.set_rotation(20)
        item.set_size(font_size)
    ax.legend_ and ax.legend_.remove()

    plt.tight_layout()

    if output_path is None:
        output_path = f"{run_name}_barplot.pdf"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved bar plot to {output_path}")
    plt.close(fig)


def barplot_lddt(
    df: pd.DataFrame,
    run_name: str,
    tools: Optional[List[str]] = None,
    lddt_col: str = "lddt",
    output_path: Optional[str] = None,
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool_name_map: Optional[dict] = None,
) -> None:
    """Create a bar plot for a single run showing LDDT bars, with
    mean runtime as a label over each bar.

    Parameters
    ----------
    df:
        Merged dataframe as returned by ``make_merged_df`` from summarize.py.
    run_name:
        Name used for the default output filename.
    tools:
        Ordered list of tool names to include. If *None* all tools in *df* are
        used, sorted by mean LDDT descending.
    output_path:
        File path to save the figure (PDF). If *None* the figure is saved as
        ``<run_name>_barplot_lddt.pdf`` in the current working directory.
    title:
        Title for the figure. Defaults to *run_name*.
    font_size:
        Base font size. Defaults to 10.
    tool_name_map:
        Mapping from dataframe tool names to display names shown on the plot.
    """
    sns.set_style("whitegrid")

    available = df["tool"].unique().tolist()

    if tools is None:
        candidates = available
    else:
        candidates = [t for t in tools if t in available]

    # Sort by mean LDDT descending
    order = (
        df[df["tool"].isin(candidates)].groupby("tool")[lddt_col].mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    if not order:
        raise ValueError("No matching tools found in the dataframe.")

    data = df[df["tool"].isin(order)].copy()
    means = data.groupby("tool")[[lddt_col, "s"]].mean()

    # Apply display-name mapping
    if tool_name_map:
        display_order = [tool_name_map.get(t, t) for t in order]
        data["tool"] = data["tool"].map(lambda t: tool_name_map.get(t, t))
        means.index = [tool_name_map.get(t, t) for t in means.index]
    else:
        display_order = order

    n_tools = len(display_order)
    font_size = font_size if font_size is not None else 10
    fig, ax = plt.subplots(figsize=(max(5, n_tools * 1.2), 4.5))
    fig.suptitle(title if title is not None else run_name, fontsize=font_size + 2, fontweight="bold")

    # Assign colors by stable alphabetical order so the same tool always
    # gets the same color regardless of the display ordering.
    stable_order = sorted(display_order)
    palette_deep = sns.color_palette("deep", n_tools)
    color_deep = {t: palette_deep[i] for i, t in enumerate(stable_order)}

    sns.barplot(
        data=data, x="tool", y=lddt_col, hue="tool", order=display_order,
        err_kws={"color": ".4", "linewidth": 1.5}, capsize=0.2,
        linewidth=1, edgecolor=".6",
        palette=color_deep, alpha=0.8, ax=ax,
    )

    # Runtime labels centred on each bar
    y_min = means[lddt_col].min()
    y_mid = means[lddt_col].max()
    label_y = y_min + (y_mid - y_min) / 2

    for tool in display_order:
        seconds = means.loc[tool, "s"]
        hours = seconds / 3600
        if hours >= 1:
            label = f"{hours:.2f} h"
        elif seconds >= 60:
            label = f"{seconds / 60:.1f} m"
        else:
            label = f"{seconds:.0f} s"
        ax.text(
            tool, label_y, label,
            ha="center", va="center",
            color="black", fontsize=font_size * 0.85,
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
        )

    ax.set_xlabel("")
    ax.set_ylabel(lddt_col, fontsize=font_size)
    ax.set_ylim(bottom=max(0, y_min - 0.05))
    for item in ax.get_xticklabels():
        item.set_rotation(20)
        item.set_size(font_size)
    ax.legend_ and ax.legend_.remove()

    plt.tight_layout()

    if output_path is None:
        output_path = f"{run_name}_barplot_{lddt_col}.pdf"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved LDDT bar plot to {output_path}")
    plt.close(fig)
