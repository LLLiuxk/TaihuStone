import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _hex_to_rgb01(hex_str: str):
    s = hex_str.strip().lstrip("#")
    if len(s) != 6:
        raise ValueError(f"Invalid hex color: {hex_str}")
    r = int(s[0:2], 16) / 255.0
    g = int(s[2:4], 16) / 255.0
    b = int(s[4:6], 16) / 255.0
    return np.array([r, g, b], dtype=float)


def _find_kt_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        name = str(c).lower()
        if "kt_weight" in name or "kt weights" in name or "kt_weights" in name:
            return c
    return df.columns[0]


def _short_name(col: str) -> str:
    s = str(col)
    if "_" in s:
        return s.split("_")[-1]
    return s


def _format_kt_label(v: object) -> str:
    s = str(v).strip()
    if not s:
        return s
    parts = [p.strip() for p in s.split(",")]
    out = []
    for p in parts:
        try:
            out.append(f"{float(p):.2f}")
        except Exception:
            out.append(p)
    return ",".join(out)


def _unified_tone_colors(n: int):
    """Sample n colors by interpolating user-provided anchor colors."""
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot translucency_final_nama.xlsx with one line per KT_weights row.")
    parser.add_argument("--input", required=True, help="xlsx path")
    parser.add_argument("--output", required=True, help="output png path")
    parser.add_argument("--title", default="Translucency Final Nama", help="plot title")
    args = parser.parse_args()

    df = pd.read_excel(args.input)
    if df.shape[1] < 2:
        raise ValueError("Need at least 2 columns: KT_weights + metrics.")

    kt_col = _find_kt_col(df)
    kt_idx = list(df.columns).index(kt_col)

    metric_cols = list(df.columns)[kt_idx + 1 :]
    if not metric_cols:
        raise ValueError("No metric columns found after KT_weights column.")

    # Keep numeric metrics only; no manual selection.
    num_df = df[metric_cols].apply(pd.to_numeric, errors="coerce")
    valid_metric_cols = [c for c in metric_cols if not num_df[c].isna().all()]
    if not valid_metric_cols:
        raise ValueError("No numeric metric columns found after KT_weights.")
    num_df = num_df[valid_metric_cols]

    # One line per row (KT_weights group label).
    labels = df[kt_col].apply(_format_kt_label)
    keep_rows = ~(num_df.isna().all(axis=1))
    num_df = num_df.loc[keep_rows].copy()
    labels = labels.loc[keep_rows]
    if num_df.empty:
        raise ValueError("No numeric rows available for plotting.")

    # Column-wise normalization with original direct factor.
    # factor = 1/max_abs (previous version behavior).
    max_abs = num_df.abs().max(axis=0)
    factors = max_abs.copy()
    factors = factors.apply(lambda v: (1.0 / v) if pd.notna(v) and abs(v) > 1e-12 else 1.0)
    norm_df = num_df.mul(factors, axis=1)

    x = np.arange(len(valid_metric_cols))
    # Fixed 16:9 canvas as requested.
    plt.figure(figsize=(16, 9), dpi=140)
    n_lines = norm_df.shape[0]
    line_colors = _unified_tone_colors(n_lines)
    for i in range(norm_df.shape[0]):
        row = norm_df.iloc[i].to_numpy(dtype=float)
        color = line_colors[i]
        plt.plot(
            x,
            row,
            marker="o",
            markersize=4,   # smaller point marker
            linewidth=1.8,
            alpha=0.95,
            color=color,
            label=labels.iloc[i]
        )

    # Use current table headers directly for x-axis labels.
    x_labels = [str(c).strip() for c in valid_metric_cols]
    plt.xticks(x, x_labels, rotation=25, ha="right")
    plt.ylabel("Normalized Value (column-wise max->1)")
    plt.xlabel("Metrics")
    plt.title(args.title)
    plt.grid(True, linestyle="--", alpha=0.35)

    # Legend on right: line labels.
    plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False, title=str(kt_col))

    # Scale notes on right side.
    note_lines = []
    for c in valid_metric_cols:
        f = factors[c]
        note_lines.append(f"{_short_name(c)} x{f:.3g}")
    note_text = "Scale factors (value * factor):\n" + "\n".join(note_lines)
    plt.gcf().text(0.82, 0.98, note_text, va="top", ha="left", fontsize=8)

    plt.tight_layout(rect=[0, 0, 0.80, 1])

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(args.output)
    print(f"[OK] Saved plot: {args.output}")
    print(f"[INFO] KT column: {kt_col}")
    print(f"[INFO] Metric columns count: {len(valid_metric_cols)}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}")
        sys.exit(1)
