"""Bar plots for benchmark results of a single run with multiple tools."""

import os
import re
import zlib
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Colors are a function of the tool name in the data only, so a tool has the
# same color in every plot, whatever other tools are shown next to it.
PALETTE = sns.color_palette("deep") + sns.color_palette("tab20b")[0::4]

# Fixed palette slots of known tools, so that they never share a color.
TOOL_COLOR_INDEX = {
    "learnMSA": 0,
    "learnMSA2": 1,
    "learnMSA_struct": 4,
    "learnMSA2_struct": 9,
    "famsa": 2,
    "clustalo": 5,
    "kalign": 6,
    "mafft": 7,
    "mafft_3Di": 8,
    "muscle_super5": 3,
    "t-coffee": 10,
    "foldmason": 11,
}


def tint(color: Tuple[float, float, float],
         factor: float = 0.45) -> Tuple[float, float, float]:
    """``color`` mixed with white, keeping ``factor`` of the color. Used
    instead of alpha so that overlaid bars stay opaque and vivid."""
    return tuple(1 - (1 - c) * factor for c in color)


def tool_color(tool: str) -> Tuple[float, float, float]:
    """The color of a tool, given by its name in the data. Unknown tools get
    a color from a stable hash of their name, which may coincide with the
    color of another tool."""
    index = TOOL_COLOR_INDEX.get(tool)
    if index is None:
        index = zlib.crc32(tool.encode()) % len(PALETTE)
    return PALETTE[index]


def tool_palettes(
    tools: List[str], tool_name_map: Optional[dict] = None,
) -> Tuple[Dict[str, tuple], Dict[str, tuple]]:
    """Tinted and full colors of ``tools`` keyed by their display names."""
    light, full = {}, {}
    for tool in tools:
        name = tool_name_map.get(tool, tool) if tool_name_map else tool
        full[name] = tool_color(tool)
        light[name] = tint(full[name])
    return light, full


def barplot(
    df: pd.DataFrame,
    run_name: str,
    tools: Optional[List[str]] = None,
    output_path: Optional[str] = None,
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool_name_map: Optional[dict] = None,
    time_col: str = "s",
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
    time_col:
        Runtime column for the labels, e.g. "s", "tool_s" or "runtime".
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
    means = data.groupby("tool")[["SP-Score", "TC", time_col]].mean()

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

    color_light, color_full = tool_palettes(order, tool_name_map)

    # SP-Score bars (faded)
    sns.barplot(
        data=data, x="tool", y="SP-Score", hue="tool", order=display_order,
        err_kws={"color": ".4", "linewidth": 1.5}, capsize=0.2,
        linewidth=1, edgecolor=".6",
        palette=color_light, ax=ax,
    )

    # TC bars (solid) overlaid
    sns.barplot(
        data=data, x="tool", y="TC", hue="tool", order=display_order,
        errorbar=None,
        palette=color_full, ax=ax,
    )

    # Runtime labels centred on each bar group
    y_min = means["TC"].min()
    y_mid = means["SP-Score"].max()
    label_y = y_min + (y_mid - y_min) / 2

    for tool in display_order:
        seconds = means.loc[tool, time_col]
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


def score_label(col: str) -> str:
    """Display name of a score column, e.g. ``lddt_muscle`` -> ``LDDT Muscle``."""
    words = []
    for word in col.split("_"):
        if "lddt" in word.lower():
            word = re.sub("lddt", "LDDT", word, flags=re.IGNORECASE)
        else:
            word = word[:1].upper() + word[1:]
        words.append(word)
    return " ".join(words)


def barplot_lddt(
    df: pd.DataFrame,
    run_name: str,
    tools: Optional[List[str]] = None,
    lddt_cols: Sequence[str] = ("lddt_muscle", "lddt"),
    output_path: Optional[str] = None,
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool_name_map: Optional[dict] = None,
    time_col: str = "s",
) -> None:
    """Create a bar plot for a single run showing LDDT bars, with
    mean runtime as a label over each bar.

    With two scores, the first (assumed to be the larger one) is drawn as
    faded bars and the second as solid bars on top, like SP-Score and TC in
    ``barplot``.

    Parameters
    ----------
    df:
        Merged dataframe as returned by ``make_merged_df`` from summarize.py.
    run_name:
        Name used for the default output filename.
    tools:
        Ordered list of tool names to include. If *None* all tools in *df* are
        used, sorted by mean of the first score descending.
    lddt_cols:
        One or two score columns to plot. The first one is assumed to be
        larger than the second.
    output_path:
        File path to save the figure (PDF). If *None* the figure is saved as
        ``<run_name>_barplot_<lddt_cols joined by _>.pdf`` in the current
        working directory.
    title:
        Title for the figure. Defaults to *run_name*.
    font_size:
        Base font size. Defaults to 10.
    tool_name_map:
        Mapping from dataframe tool names to display names shown on the plot.
    time_col:
        Runtime column for the labels, e.g. "s", "tool_s" or "runtime".
    """
    sns.set_style("whitegrid")

    lddt_cols = [lddt_cols] if isinstance(lddt_cols, str) else list(lddt_cols)
    if not 1 <= len(lddt_cols) <= 2:
        raise ValueError("lddt_cols must contain one or two score columns.")
    top_col, bottom_col = lddt_cols[0], lddt_cols[-1]

    available = df["tool"].unique().tolist()

    if tools is None:
        candidates = available
    else:
        candidates = [t for t in tools if t in available]

    # Sort by the mean of the first score descending
    order = (
        df[df["tool"].isin(candidates)].groupby("tool")[top_col].mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    if not order:
        raise ValueError("No matching tools found in the dataframe.")

    data = df[df["tool"].isin(order)].copy()
    means = data.groupby("tool")[lddt_cols + [time_col]].mean()

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

    color_light, color_full = tool_palettes(order, tool_name_map)

    two_scores = len(lddt_cols) == 2

    # First (larger) score; faded if a second one is drawn on top
    sns.barplot(
        data=data, x="tool", y=top_col, hue="tool", order=display_order,
        err_kws={"color": ".4", "linewidth": 1.5}, capsize=0.2,
        linewidth=1, edgecolor=".6",
        palette=color_light if two_scores else color_full, ax=ax,
    )

    # Second (smaller) score, solid, overlaid
    if two_scores:
        sns.barplot(
            data=data, x="tool", y=bottom_col, hue="tool", order=display_order,
            errorbar=None,
            palette=color_full, ax=ax,
        )

    # Runtime labels centred on each bar
    y_min = means[bottom_col].min()
    y_mid = means[top_col].max()
    label_y = y_min + (y_mid - y_min) / 2

    for tool in display_order:
        seconds = means.loc[tool, time_col]
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
    if two_scores:
        ylabel = (f"{score_label(top_col)} (fade) / "
                  f"{score_label(bottom_col)} (solid)")
    else:
        ylabel = score_label(top_col)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.set_ylim(bottom=max(0, y_min - 0.05))
    for item in ax.get_xticklabels():
        item.set_rotation(20)
        item.set_size(font_size)
    ax.legend_ and ax.legend_.remove()

    plt.tight_layout()

    if output_path is None:
        output_path = f"{run_name}_barplot_{'_'.join(lddt_cols)}.pdf"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved LDDT bar plot to {output_path}")
    plt.close(fig)
