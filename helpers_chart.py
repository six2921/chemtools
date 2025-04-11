"""
Reusable plotting helpers for Jupyter.

• plot_histogram()   : histogram with optional hue‑split + count labels
• chart_num()        : interactive histogram UI  (numeric columns)
• plot_stacked_bar() : simple proportion bar chart
• chart_category()   : interactive stacked‑bar UI (≤10 unique values)
"""

from __future__ import annotations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from ipywidgets import (
    Dropdown, FloatText, Button, VBox, HBox, Layout, Output
)
from IPython.display import display, clear_output


# ────────────────────────────────────────────────────────────────
# Histogram (single function)
# ────────────────────────────────────────────────────────────────
def plot_histogram(
    df: pd.DataFrame,
    x_col: str,
    hue_col: str,
    vmin: float,
    vmax: float,
    bin_width: float,
    ratio: float = 1.0,
    save_path: str | None = None,
) -> None:
    """Draw a histogram; annotate counts; optional PNG export."""
    data = df[(df[x_col] >= vmin) & (df[x_col] <= vmax)]
    fig, ax = plt.subplots(figsize=(6 * ratio, 6))
    bins = np.arange(vmin, vmax + bin_width, bin_width)

    if hue_col != "None":
        sns.histplot(
            data,
            x=x_col,
            hue=hue_col,
            bins=bins,
            multiple="dodge",
            shrink=0.8,
            palette="tab10",
            ax=ax,
        )
    else:
        sns.histplot(data[x_col], bins=bins, color="tab:blue", shrink=0.8, ax=ax)

    for p in ax.patches:
        if (h := p.get_height()) > 0:
            ax.text(
                p.get_x() + p.get_width() / 2,
                h,
                int(h),
                ha="center",
                va="bottom",
            )

    ax.set_xticks(bins)
    ax.set_xlabel(x_col, fontsize=14, fontweight="bold")
    ax.set_ylabel("Count")
    ax.yaxis.grid(True, ls="--", lw=0.5)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()


# ────────────────────────────────────────────────────────────────
# Interactive histogram UI (numeric)
# ────────────────────────────────────────────────────────────────
def chart_num(df: pd.DataFrame) -> None:
    """Launch interactive histogram builder for numeric columns."""
    num_cols = df.select_dtypes(["number"]).columns.tolist()
    hue_opts = ["None"] + [c for c in df.columns if df[c].nunique() <= 10]

    w_x      = Dropdown(options=num_cols, value=num_cols[0], description="Column:", layout=Layout(width="70%"))
    w_hue    = Dropdown(options=hue_opts, value="None",     description="Hue:",    layout=Layout(width="70%"))
    w_min    = FloatText(description="Min:",       layout=Layout(width="70%"))
    w_max    = FloatText(description="Max:",       layout=Layout(width="70%"))
    w_bin    = FloatText(value=1.0, description="Bin width:", layout=Layout(width="70%"))
    w_ratio  = FloatText(value=1.0, description="Image ratio:", layout=Layout(width="70%"))
    btn_save = Button(description="Export PNG", layout=Layout(width="70%"))
    out      = Output()

    def _update_range(_=None):
        col = w_x.value
        w_min.value, w_max.value = float(df[col].min()), float(df[col].max())
    w_x.observe(_update_range, names="value")
    _update_range()

    def _draw(_=None):
        with out:
            clear_output(wait=True)
            plot_histogram(
                df, w_x.value, w_hue.value,
                w_min.value, w_max.value,
                w_bin.value, w_ratio.value
            )
    for w in (w_x, w_hue, w_min, w_max, w_bin, w_ratio):
        w.observe(_draw, names="value")

    def _save(_):
        path = f"{w_x.value}.png"
        plot_histogram(
            df, w_x.value, w_hue.value,
            w_min.value, w_max.value,
            w_bin.value, w_ratio.value,
            save_path=path
        )
        print(f"✅  Saved to {path}")
    btn_save.on_click(_save)

    controls = VBox([w_x, w_hue, w_min, w_max, w_bin, w_ratio, btn_save],
                    layout=Layout(align_items="flex-start"))
    display(HBox([controls, out], layout=Layout(align_items="center")))


# ────────────────────────────────────────────────────────────────
# Proportion bar chart (single function)
# ────────────────────────────────────────────────────────────────
def plot_stacked_bar(df: pd.DataFrame, col: str):
    """Return a figure with proportion bar chart and count labels."""
    counts = df[col].astype(str).value_counts()
    prop   = counts / counts.sum()

    fig, ax = plt.subplots(figsize=(3, 3))
    bars = ax.bar(prop.index, prop.values, color=plt.cm.tab10(range(len(prop))))

    for bar, p, n in zip(bars, prop.values, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, p / 2,
                f"{p*100:.1f}%\n({n})", ha="center", va="center")

    ax.set_ylabel("Proportion")
    ax.set_ylim(0, 1)
    ax.set_xlabel(col, fontsize=14, fontweight="bold")
    ax.yaxis.grid(True, ls="--", lw=0.5)
    plt.show()
    return fig


# ────────────────────────────────────────────────────────────────
# Interactive stacked‑bar UI (categorical ≤10 unique)
# ────────────────────────────────────────────────────────────────
def chart_category(df: pd.DataFrame):
    """Launch interactive stacked‑bar chart builder for low‑cardinality columns."""
    cat_cols = [c for c in df.columns if df[c].nunique() <= 10]
    if not cat_cols:
        print("⚠️  No columns with ≤10 unique values.")
        return

    w_col    = Dropdown(options=cat_cols, description="Column:", layout=Layout(width="70%"))
    btn_save = Button(description="Export PNG", layout=Layout(width="70%"))
    out      = Output()

    def _draw(_=None):
        with out:
            clear_output(wait=True)
            fig = plot_stacked_bar(df, w_col.value)
            btn_save.on_click(lambda *_: fig.savefig(f"{w_col.value}.png", dpi=300, bbox_inches="tight"))

    w_col.observe(_draw, names="value")
    _draw()

    display(HBox([VBox([w_col, btn_save]), out], layout=Layout(align_items="center")))

