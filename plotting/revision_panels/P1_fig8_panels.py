"""Panels added to or replacing parts of Fig. 8 during revision.

Style comes from fig8_style.py (Arial 8 pt, no titles, NPG palette, PDF and
PNG at 300 dpi). The Venn follows 04_venn.py and the forest plots follow
09_cox_forest.py, so the new panels sit beside the published ones without
restyling.

Outputs, each as .pdf and .png:
  8A_alt__venn_n1500_strict     three-way overlap at 1,500 sender cells with
                                the stricter detection filter (Fig. S13)
  8_new__sec61g_cn_adjusted     SEC61G survival with copy number added
                                stepwise (Fig. 8O)
  8_new__locus_specificity      each 7p11.2 gene substituted into the same
                                model (Fig. 8P)
  8_new__spatial_depth_section  spatial endpoint after per-section adjustment
                                for sequencing depth (Fig. 8N)

Usage: python P1_fig8_panels.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fig8"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn3, venn3_circles

from fig8_style import (setup, COL, save_panel, sig_stars, add_subtitle,
                        style_axes)
setup()

# Journal convention: italic P, scientific notation. matplotlib's mathtext
# defaults to DejaVu, which would render the italic P in a different face from
# the Arial used everywhere else, so point mathtext at Arial explicitly.
plt.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
})


def psci(p, thresh=0.01):
    """Italic P; scientific notation below `thresh`, three decimals above."""
    if p is None or not np.isfinite(p):
        return r"$\it{P}$ = n.a."
    if p <= 0:
        return r"$\it{P}$ < 10$^{-308}$"
    if p < thresh:
        e = int(np.floor(np.log10(p)))
        return f"$\\it{{P}}$ = {p/10**e:.1f}\u00d710$^{{{e}}}$"
    return f"$\\it{{P}}$ = {p:.3f}"

RES = Path("${PROJECT_ROOT}/results")
OUT = Path("${WORK_ROOT}/revision_panels")
OUT.mkdir(parents=True, exist_ok=True)

# Alternative to Fig. 8A: 1,500 sender cells, stricter detection filter.
# Identical to 04_venn.py except for the input directory and the filter:
# n_det >= 5 becomes detection in >= 20% of cells, the point at which the
# correlation between the goal-shift score and detection frequency first
# becomes indistinguishable from zero.
TRANSITIONS = {
    "macro_spp1_to_c1qc": ("Macro_SPP1→C1QC", COL["venn_macro"],
                           RES / "fig8_geneformer/perturb_rent_macro", 1500),
    "mal_mp3_to_mp1": ("Mal_MP3→MP1", COL["venn_mal"],
                       RES / "fig8_geneformer/perturb_rent_mal", 1500),
    "neu_osm_priming_to_low": ("Neu_OSM_priming→low", COL["venn_neu"],
                               RES / "fig8_geneformer/perturb_neuFULL_n1496",
                               1496),
}
FDR, FRAC = 0.05, 0.20


def load_filtered(t):
    lab, col, root, ncell = TRANSITIONS[t]
    d = pd.read_csv(root / t / f"{t}_stats.csv", index_col=0)
    d = d[(d["Sig"] == 1) & (d["Shift_to_goal_end"] > 0)
          & (d["Goal_end_FDR"] < FDR)
          & (d["N_Detections"] >= int(FRAC * ncell))].copy()
    return d.sort_values("Shift_to_goal_end", ascending=False
                         ).reset_index(drop=True)


filt = {t: load_filtered(t) for t in TRANSITIONS}
for t, d in filt.items():
    print(f"  {t}: {len(d)} candidates pass filter")


def gset(d, n):
    return set(d.head(n)["Gene_name"])


def make_venn_with_bar(top_n, out_stem):
    a, b, c = (gset(filt[t], top_n) for t in TRANSITIONS)
    only_a, only_b, only_c = a - b - c, b - a - c, c - a - b
    ab, ac, bc = (a & b) - c, (a & c) - b, (b & c) - a
    abc = a & b & c

    fig = plt.figure(figsize=(5.6, 3.4))
    gs = GridSpec(1, 2, width_ratios=[3.5, 1.0], wspace=0.05, figure=fig)
    ax_v = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])

    subsets = (len(only_a), len(only_b), len(ab), len(only_c),
               len(ac), len(bc), len(abc))
    labels = [TRANSITIONS[t][0] for t in TRANSITIONS]
    cols_v = tuple(TRANSITIONS[t][1] for t in TRANSITIONS)
    v = venn3(subsets=subsets, set_labels=labels, set_colors=cols_v,
              alpha=0.55, ax=ax_v)
    venn3_circles(subsets=subsets, linewidth=0.6, ax=ax_v, color="black")
    for lbl in list(v.set_labels) + list(v.subset_labels):
        if lbl is not None:
            lbl.set_fontsize(7)
    ax_v.text(0.5, -0.07, "3-way: none" if not abc else
              f"3-way: {', '.join(sorted(abc))}",
              ha="center", fontsize=6, transform=ax_v.transAxes)
    add_subtitle(ax_v, f"top-{top_n} per transition · 1,500 cells · "
                       f"detected in ≥{FRAC:.0%} of cells · FDR < {FDR}")

    bar_short = ["Macro", "Mal", "Neu"]
    bar_vals = [len(filt[t]) for t in TRANSITIONS]
    bar_cols = [TRANSITIONS[t][1] for t in TRANSITIONS]
    bars = ax_b.barh(bar_short, bar_vals, color=bar_cols, edgecolor="black",
                     linewidth=0.4)
    for bar, val in zip(bars, bar_vals):
        ax_b.text(val + max(bar_vals) * 0.02,
                  bar.get_y() + bar.get_height() / 2, str(val), va="center",
                  fontsize=6)
    style_axes(ax_b, xlabel="# genes passing filter", ylabel="")
    ax_b.invert_yaxis()

    fig.tight_layout()
    save_panel(fig, str(out_stem))
    plt.close(fig)
    pool = (a & b) | (a & c) | (b & c)
    return abc, pool


abc, pool = make_venn_with_bar(200, OUT / "8A_alt__venn_n1500_strict")
print(f"  ≥2-hit pool: {len(pool)} genes; 3-way: {len(abc)}")

# Two forest plots, drawn as in 09_cox_forest.py
def forest(df, gene_col, label_fmt, xlabel, subtitle, out_stem,
           highlight=None, figsize=(4.6, 2.4), xticks=(0.9, 1.0, 1.2, 1.5, 2.0)):
    fig, ax = plt.subplots(figsize=figsize)
    for i, row in df.iterrows():
        base = COL["hr_up"] if row["HR"] >= 1 else COL["hr_down"]
        color = base if (highlight is None or row[gene_col] == highlight) \
            else "#B0B0B0"
        ax.hlines(i, row["CI_low"], row["CI_high"], color=color, lw=0.8,
                  alpha=0.9)
        ax.plot(row["HR"], i, marker="s", color=color, markersize=4.0,
                markeredgecolor="white", markeredgewidth=0.4)
    ax.axvline(1.0, color="black", lw=0.4, ls="--", alpha=0.7)
    ax.set_xscale("log")
    xmin = min(df["CI_low"].min(), 0.9) * 0.94
    data_xmax = max(df["CI_high"].max(), 1.5)
    ax.set_xlim(xmin, data_xmax * 3.2)
    for i, row in df.iterrows():
        ax.text(data_xmax * 1.06, i, label_fmt(row), ha="left", va="center",
                fontsize=6, color="black")
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df[gene_col],
                       style="italic" if gene_col == "gene" else "normal")
    ax.invert_yaxis()
    ax.set_xticks(list(xticks))
    ax.set_xticklabels([str(t) for t in xticks])
    ax.minorticks_off()
    style_axes(ax, xlabel=xlabel, ylabel="")
    add_subtitle(ax, subtitle)
    fig.tight_layout()
    save_panel(fig, str(out_stem))
    plt.close(fig)


CX = pd.read_csv(RES / "cn_confounding/sec61g_cox_ladder.csv")
forest(CX, "model",
       lambda r: f"HR = {r.HR:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f})  "
                 f"{psci(r.p)} {sig_stars(r.p)}",
       "Hazard ratio (SEC61G expression, per 1 SD)",
       f"TCGA-LUAD n={int(CX.n.iloc[0])} ({int(CX.events.iloc[0])} events) · "
       "covariates added stepwise",
       OUT / "8_new__sec61g_cn_adjusted", figsize=(5.0, 2.2))

LS = pd.read_csv(RES / "cn_confounding/locus_specificity.csv"
                 ).sort_values("HR").reset_index(drop=True)
forest(LS, "gene",
       lambda r: f"HR = {r.HR:.2f} ({r.CI_low:.2f}–{r.CI_high:.2f})  "
                 f"{psci(r.p)} {sig_stars(r.p)}",
       "Hazard ratio (per 1 SD, identical adjusted model)",
       "each 7p11.2 gene substituted into: gene + 7p11.2 CN + age + stage + sex",
       OUT / "8_new__locus_specificity", highlight="SEC61G",
       figsize=(5.0, 2.4))

# Spatial endpoint: per-section adjustment for sequencing depth
S = {t: pd.read_csv(RES / f"spatial_depth_per_section/{t}_SEC61G_by_section.csv")
     for t in ("discovery", "validation")}
Q = {t: pd.read_csv(RES / f"spatial_depth_per_section/{t}_random.csv")
     for t in ("discovery", "validation")}

fig, ax = plt.subplots(figsize=(3.6, 2.0))
rng = np.random.default_rng(0)
ypos = {"discovery": 1, "validation": 0}
for t, y0 in ypos.items():
    lo, hi = np.percentile(Q[t].beta_med, [5, 95])
    ax.fill_betweenx([y0 - 0.26, y0 + 0.26], lo, hi, color="#CCCCCC",
                     alpha=0.45, lw=0)
    b = S[t].beta.values
    ax.plot(b, y0 + rng.uniform(-0.10, 0.10, len(b)), "o", ms=3.6,
            color=COL["venn_macro"], mec="white", mew=0.4, zorder=3)
    ax.plot([np.median(b)], [y0], "|", ms=13, mew=1.6, color=COL["hr_up"],
            zorder=4)
ax.axvline(0, color="black", lw=0.4, ls="--", alpha=0.7)
ax.set_yticks([1, 0])
ax.set_yticklabels([f"Discovery\n({len(S['discovery'])} sections)",
                    f"Validation\n({len(S['validation'])} sections)"])
ax.set_ylim(-0.6, 1.6)
style_axes(ax, xlabel="SEC61G effect per section, adjusted for spot depth",
           ylabel="")
add_subtitle(ax, "grey band = 5–95% of expression-matched random genes")
fig.tight_layout()
save_panel(fig, str(OUT / "8_new__spatial_depth_section"))
plt.close(fig)

print("\nwritten to", OUT)
