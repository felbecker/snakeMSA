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
    "muscle_3d": 12,
    "muscle_3d_super7": 13,
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


def replica_means(data: pd.DataFrame, col: str) -> Dict[str, np.ndarray]:
    """Per tool, the mean of ``col`` over all rows (families) of each replica.

    Only tools with at least two replicas are returned, so a tool without
    replicas, or a dataframe without a ``replica`` column, gets no error bar.
    """
    if "replica" not in data.columns:
        return {}
    per_replica = data.groupby(["tool", "replica"])[col].mean().dropna()
    result = {}
    for tool, values in per_replica.groupby(level="tool"):
        if len(values) >= 2:
            result[tool] = values.to_numpy()
    return result


def draw_replica_errors(
    ax: plt.Axes, data: pd.DataFrame, col: str, display_order: List[str],
) -> None:
    """Draw, at the top of each bar of ``col``, the error of the replicas: the
    sample standard deviation of the per-replica means over all families, as
    an error bar, together with the individual replica means as dots."""
    bar_heights = data.groupby("tool")[col].mean()
    per_tool = replica_means(data, col)
    for i, tool in enumerate(display_order):
        values = per_tool.get(tool)
        if values is None:
            continue
        ax.errorbar(
            i, bar_heights[tool], yerr=np.std(values, ddof=1), fmt="none",
            ecolor=".2", elinewidth=1.5, capsize=4, capthick=1.5, zorder=4,
        )
        offsets = np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 \
            else [0]
        ax.scatter(
            i + np.asarray(offsets), values, s=14, color=".2",
            edgecolor="white", linewidth=0.5, zorder=5,
        )


def sorted_tools(
    df: pd.DataFrame, candidates: List[str], col: str, sort_by_score: bool,
) -> List[str]:
    """``candidates`` sorted by the mean of ``col`` descending, or in their
    given order if not ``sort_by_score``."""
    if not sort_by_score:
        return list(candidates)
    return (
        df[df["tool"].isin(candidates)].groupby("tool")[col].mean()
        .sort_values(ascending=False)
        .index.tolist()
    )


def format_runtime(seconds: float) -> str:
    """A runtime label in the largest fitting unit (s, m or h)."""
    if seconds >= 3600:
        return f"{seconds / 3600:.2f} h"
    if seconds >= 60:
        return f"{seconds / 60:.1f} m"
    return f"{seconds:.0f} s"


def draw_runtime(
    ax: plt.Axes,
    seconds: pd.Series,
    display_order: List[str],
    label_y: float,
    font_size: float,
    runtime_line: bool = False,
) -> None:
    """Show the mean runtime of each tool.

    By default, as a text label at height ``label_y`` over each bar. With
    ``runtime_line``, as red dots connected by a red line on a second y axis
    on the right hand side, in the unit that fits the largest runtime.
    """
    seconds = seconds.loc[display_order]
    if not runtime_line:
        for tool in display_order:
            ax.text(
                tool, label_y, format_runtime(seconds[tool]),
                ha="center", va="center",
                color="black", fontsize=font_size * 0.85,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            )
        return

    longest = seconds.max()
    if longest >= 3600:
        unit, scale = "h", 3600
    elif longest >= 60:
        unit, scale = "min", 60
    else:
        unit, scale = "s", 1
    rt_ax = ax.twinx()
    # The grid of the score axis is enough; a second one would not align
    rt_ax.grid(False)
    rt_ax.plot(
        np.arange(len(display_order)), seconds.to_numpy() / scale,
        color="red", marker="o", markersize=6, linewidth=1.5, zorder=6,
    )
    rt_ax.set_ylim(bottom=0)
    rt_ax.set_ylabel(f"Runtime ({unit})", color="red", fontsize=font_size)
    rt_ax.tick_params(axis="y", colors="red", labelsize=font_size)
    rt_ax.spines["right"].set_color("red")


def barplot(
    df: pd.DataFrame,
    run_name: str,
    tools: Optional[List[str]] = None,
    output_path: Optional[str] = None,
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool_name_map: Optional[dict] = None,
    time_col: str = "s",
    runtime_line: bool = False,
    sort_by_score: bool = True,
) -> None:
    """Create a bar plot for a single run overlaying SP-Score (faded) and TC
    (solid) bars, with mean runtime as a label over each bar group.

    Parameters
    ----------
    df:
        Merged dataframe as returned by ``make_merged_df`` from summarize.py,
        already filtered to this run.
        With a ``replica`` column (``--average-replicas``), error bars show
        the standard deviation of the per-replica means over all families.
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
    runtime_line:
        Instead of runtime labels, plot the runtimes as red dots connected
        by a red line on a second y axis on the right hand side.
    sort_by_score:
        Sort the tools by score. If *False*, plot them in the order of
        *tools* (all tools in *df* in order of appearance if *None*).
    """
    sns.set_style("whitegrid")

    available = df["tool"].unique().tolist()

    if tools is None:
        candidates = available
    else:
        candidates = [t for t in tools if t in available]

    order = sorted_tools(df, candidates, "TC", sort_by_score)

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
        errorbar=None,
        linewidth=1, edgecolor=".6",
        palette=color_light, ax=ax,
    )

    # TC bars (solid) overlaid
    sns.barplot(
        data=data, x="tool", y="TC", hue="tool", order=display_order,
        errorbar=None,
        palette=color_full, ax=ax,
    )

    draw_replica_errors(ax, data, "SP-Score", display_order)
    draw_replica_errors(ax, data, "TC", display_order)

    # Runtime labels centred on each bar group
    y_min = means["TC"].min()
    y_mid = means["SP-Score"].max()
    label_y = y_min + (y_mid - y_min) / 2
    draw_runtime(ax, means[time_col], display_order, label_y, font_size,
                 runtime_line)

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
    runtime_line: bool = False,
    sort_by_score: bool = True,
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
        With a ``replica`` column (``--average-replicas``), error bars show
        the standard deviation of the per-replica means over all families.
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
    runtime_line:
        Instead of runtime labels, plot the runtimes as red dots connected
        by a red line on a second y axis on the right hand side.
    sort_by_score:
        Sort the tools by score. If *False*, plot them in the order of
        *tools* (all tools in *df* in order of appearance if *None*).
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

    order = sorted_tools(df, candidates, top_col, sort_by_score)

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
        errorbar=None,
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

    for col in dict.fromkeys(lddt_cols):
        draw_replica_errors(ax, data, col, display_order)

    # Runtime labels centred on each bar
    y_min = means[bottom_col].min()
    y_mid = means[top_col].max()
    label_y = y_min + (y_mid - y_min) / 2
    draw_runtime(ax, means[time_col], display_order, label_y, font_size,
                 runtime_line)

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


#: Marker shapes of the tools drawn as single points next to the curve.
OTHER_TOOL_MARKERS = ["s", "^", "D", "v", "P", "X", "*", "<", ">"]


def draw_replica_spread(
    ax: plt.Axes, x: float, mean: float, values: np.ndarray,
) -> None:
    """Error of the replicas as in the bar plots: the standard deviation of
    the per-replica means as an error bar at ``mean``, with the replica
    means as dots."""
    if len(values) < 2:
        return
    ax.errorbar(x, mean, yerr=np.std(values, ddof=1), fmt="none",
                ecolor=".2", elinewidth=1.5, capsize=4, capthick=1.5,
                zorder=4)
    ax.scatter(np.full(len(values), x), values, s=14, color=".2",
               edgecolor="white", linewidth=0.5, zorder=5)


def missing_struct_curve(
    df: pd.DataFrame,
    score_col: str = "lddt_muscle",
    output_path: str = "missing_struct_curve.pdf",
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool: str = "learnMSA_struct",
    others: Optional[pd.DataFrame] = None,
    tool_name_map: Optional[dict] = None,
    curve_label: Optional[str] = None,
) -> None:
    """Line plot of a score over the fraction p of available structures,
    from p = 1 (all structures) on the left to p = 0 on the right.

    Parameters
    ----------
    df:
        One row per alignment with the columns ``p``, ``replica`` and
        ``score_col``.
    score_col:
        Score on the y axis.
    output_path:
        File path to save the figure (PDF).
    title:
        Title for the figure. No title if *None*.
    font_size:
        Base font size. Defaults to 10.
    tool:
        Tool of the curve, for its color and legend label.
    others:
        Other tools, drawn as one point each (own marker, in their usual
        color) at the p given in their rows, with a legend. One row per
        alignment with the columns ``tool``, ``p``, ``replica`` and
        ``score_col``.
    tool_name_map:
        Mapping from tool names to legend labels.
    curve_label:
        Legend label of the curve, instead of the mapped name of *tool*.
    """
    sns.set_style("whitegrid")
    font_size = font_size if font_size is not None else 10
    names = tool_name_map or {}

    def label(name: str) -> str:
        return names.get(name, name).replace("\n", " ")

    means = df.groupby("p")[score_col].mean().sort_index(ascending=False)
    per_replica = df.groupby(["p", "replica"])[score_col].mean()

    with_others = others is not None and len(others) > 0
    # Room for the legend next to the axes
    fig, ax = plt.subplots(figsize=(7.5 if with_others else 5, 4))
    if title is not None:
        fig.suptitle(title, fontsize=font_size + 2, fontweight="bold")
    ax.plot(means.index, means.values, color=tool_color(tool), marker="o",
            lw=2, zorder=3,
            label=curve_label if curve_label is not None else label(tool))
    for p, values in per_replica.groupby(level="p"):
        draw_replica_spread(ax, p, means[p], values.to_numpy())

    if with_others:
        order = list(dict.fromkeys(others["tool"]))
        for k, other in enumerate(order):
            rows = others[others["tool"] == other]
            mean = rows[score_col].mean()
            x = rows["p"].iloc[0]
            ax.scatter(x, mean, s=60, color=tool_color(other),
                       marker=OTHER_TOOL_MARKERS[k % len(OTHER_TOOL_MARKERS)],
                       edgecolor=".2", linewidth=0.6, zorder=6,
                       label=label(other))
            replicas = rows.groupby("replica")[score_col].mean().to_numpy()
            draw_replica_spread(ax, x, mean, replicas)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                  fontsize=font_size * 0.85, frameon=False)
    # Start with all structures (p = 1) on the left
    ax.set_xlim(1.05, -0.05)

    ax.set_xticks(means.index)
    ax.set_xticklabels([f"{p:g}" for p in means.index], fontsize=font_size)
    ax.set_xlabel("Fraction of available structures p", fontsize=font_size)
    ax.set_ylabel(score_label(score_col), fontsize=font_size)
    ax.tick_params(axis="y", labelsize=font_size)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved missing structure curve to {output_path}")
    plt.close(fig)


def format_depth(depth: float) -> str:
    """Tick label of an alignment depth, e.g. 50000 -> ``50k``."""
    for factor, suffix in ((1_000_000, "M"), (1000, "k")):
        if depth >= factor and depth % factor == 0:
            return f"{depth // factor:g}{suffix}"
    return f"{depth:g}"


def depth_curve(
    df: pd.DataFrame,
    score_col: str = "TC",
    output_path: str = "depth_curve.pdf",
    tools: Optional[Sequence[str]] = None,
    title: Optional[str] = None,
    font_size: Optional[int] = None,
    tool_name_map: Optional[dict] = None,
    logy: bool = False,
) -> None:
    """Line plot of a score over the alignment depth (number of aligned
    sequences, log scale), one curve per tool in its usual color.

    Parameters
    ----------
    df:
        One row per alignment with the columns ``tool``, ``depth``,
        ``replica`` and ``score_col``.
    score_col:
        Score on the y axis, averaged per depth over all rows (families,
        subsamples and replicas).
    output_path:
        File path to save the figure (PDF).
    tools:
        Tools to draw, in legend order. Default: all tools in *df*.
    title:
        Title for the figure. No title if *None*.
    font_size:
        Base font size. Defaults to 10.
    tool_name_map:
        Mapping from tool names to legend labels.
    logy:
        Logarithmic y axis, e.g. for runtimes.
    """
    sns.set_style("whitegrid")
    font_size = font_size if font_size is not None else 10
    names = tool_name_map or {}
    if tools is None:
        tools = list(dict.fromkeys(df["tool"]))

    fig, ax = plt.subplots(figsize=(8, 4))
    if title is not None:
        fig.suptitle(title, fontsize=font_size + 2, fontweight="bold")
    # Shift the error bars of the tools apart a little on the log x axis
    shifts = 1.04 ** (np.arange(len(tools)) - (len(tools) - 1) / 2)
    for k, tool in enumerate(tools):
        rows = df[df["tool"] == tool]
        means = rows.groupby("depth")[score_col].mean().sort_index()
        ax.plot(means.index * shifts[k], means.values,
                color=tool_color(tool),
                marker=OTHER_TOOL_MARKERS[k % len(OTHER_TOOL_MARKERS)],
                markersize=7, markeredgecolor=".2", markeredgewidth=0.6,
                lw=2, zorder=3,
                label=names.get(tool, tool).replace("\n", " "))
        per_replica = rows.groupby(["depth", "replica"])[score_col].mean()
        for depth, values in per_replica.groupby(level="depth"):
            draw_replica_spread(ax, depth * shifts[k], means[depth],
                                values.dropna().to_numpy())

    depths = np.sort(df["depth"].unique())
    ax.set_xscale("log")
    ax.set_xticks(depths)
    ax.set_xticklabels([format_depth(d) for d in depths], fontsize=font_size)
    ax.minorticks_off()
    ax.set_xlim(depths[0] / 1.4, depths[-1] * 1.4)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("Alignment depth (number of sequences)", fontsize=font_size)
    ax.set_ylabel("Runtime (s)" if score_col == "tool_s"
                  else score_label(score_col), fontsize=font_size)
    ax.tick_params(axis="y", labelsize=font_size)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=font_size * 0.85, frameon=False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved depth curve to {output_path}")
    plt.close(fig)
