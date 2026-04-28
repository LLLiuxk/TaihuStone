import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 14
})

def _hex_to_rgb01(hex_str: str):
    s = hex_str.strip().lstrip("#")
    if len(s) != 6:
        raise ValueError(f"Invalid hex color: {hex_str}")
    r = int(s[0:2], 16) / 255.0
    g = int(s[2:4], 16) / 255.0
    b = int(s[4:6], 16) / 255.0
    return np.array([r, g, b], dtype=float)

def _unified_tone_colors(n: int):
    if n <= 0:
        return []
    anchors_hex = ["FE9C94", "F6E1B6", "B6DA04", "A8D2B0", "00B4CC", "C1D8E9", "92B1D9", "B5B8E4", "D8C0EE", "C4C4C4"]
    anchors = np.array([_hex_to_rgb01(h) for h in anchors_hex], dtype=float)
    if n == 1:
        c = anchors[0]
        return [(float(c[0]), float(c[1]), float(c[2]))]

    t_anchor = np.linspace(0.0, 1.0, anchors.shape[0])
    t_query = np.linspace(0.0, 1.0, n)
    r = np.interp(t_query, t_anchor, anchors[:, 0])
    g = np.interp(t_query, t_anchor, anchors[:, 1])
    b = np.interp(t_query, t_anchor, anchors[:, 2])
    return [(float(r[i]), float(g[i]), float(b[i])) for i in range(n)]

xlsx_path = Path("translucency_final_nama.xlsx")
sheet_name = "translucency_summary_merged"
df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

original_labels = [
    "Only Dir.",
    "Only Loc.",
    "Only Len.",
    "w/o Ang.",
    "Only Ang.",
    "Average",
    "Default",
    "w/o Len.",
    "w/o Loc.",
    "w/o Dir.",
]

# Required x-axis label order
requested_order = [
    "Default",
    "Average",
    "Only Ang.",
    "Only Len.",
    "Only Loc.",
    "Only Dir.",
    "w/o Ang.",
    "w/o Len.",
    "w/o Loc.",
    "w/o Dir."
]

# Columns from excel
core_metrics_map = {
    "VP": "Structure VP",
    "Porosity": "Porosity",
    "Unsupported Voxel": "Illegal voxel num",
}

path_metrics_map = {
    "Straightness": "Avg paths straightness",
    "Length": "Avg paths length",
    "Inter-ratio": "Inner kernel ratio",
    "Horizontal": "Paths horiz score",
}

def scale_and_order_df(mapping_dict):
    sub_df = df[list(mapping_dict.values())].copy()
    sub_df.columns = list(mapping_dict.keys())
    sub_df.index = requested_order

    
    # Max scaling
    norm_df = sub_df.copy()
    for col in norm_df.columns:
        max_val = norm_df[col].abs().max()
        if max_val > 1e-12:
            norm_df[col] = norm_df[col] / max_val
            
    return norm_df


def plot_grouped_bar_on_ax(ax, plot_df, title, is_bottom=False, custom_colors_hex=None):
    x = np.arange(len(plot_df.index))
    n_metrics = len(plot_df.columns)
    bar_width = 0.78 / n_metrics
    
    if custom_colors_hex:
        bar_colors = [tuple(_hex_to_rgb01(h)) for h in custom_colors_hex]
    else:
        bar_colors = _unified_tone_colors(n_metrics)

    for i, metric in enumerate(plot_df.columns):
        offset = (i - (n_metrics - 1) / 2) * bar_width
        ax.bar(x + offset, plot_df[metric].values, width=bar_width, label=metric, color=bar_colors[i])

    ax.set_ylabel("Normalized value")
    if is_bottom:
        ax.set_xlabel("Weight configuration")
    else:
        ax.set_xlabel("")
        
    ax.set_title(title, pad=45)
    ax.set_xticks(x)
    
    # Both top and bottom plot will show the x labels
    ax.set_xticklabels(plot_df.index, rotation=30, ha="right")
        
    ax.set_ylim(0, 1.08)
    ax.legend(ncol=n_metrics, loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.35)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

df1 = scale_and_order_df(core_metrics_map)
df2 = scale_and_order_df(path_metrics_map)

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(11.5, 11.0))

plot_grouped_bar_on_ax(axes[0], df1, "Sensitivity Analysis of VP Weights (Core Metrics)", is_bottom=False, custom_colors_hex=["00B4CC", "FE9C94", "B6DA04"])
plot_grouped_bar_on_ax(axes[1], df2, "Sensitivity Analysis of VP Weights (Path Metrics)", is_bottom=True)

fig.tight_layout()
fig.subplots_adjust(hspace=0.70)

fig.savefig("vp_weight_sensitivity_combined.png", dpi=300, bbox_inches="tight")
fig.savefig("vp_weight_sensitivity_combined.pdf", bbox_inches="tight")
plt.close(fig)

print("[OK] Saved combined grouped bar charts successfully.")
