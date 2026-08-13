"""Fig 8A: 3-way Venn + per-transition geneset size barplot."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn3, venn3_circles

import sys
from pathlib import Path
# fig8_style lives with the figure it styles; make it importable when this
# script is run from its own directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "plotting" / "fig8"))
from fig8_style import setup, COL, save_panel, add_subtitle, style_axes
setup()

PERT = Path("${PROJECT_ROOT}/results/fig8_geneformer/perturb")
OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn")
OUT.mkdir(parents=True, exist_ok=True)

TRANSITIONS = {
    "macro_spp1_to_c1qc":     ("Macro_SPP1→C1QC", COL["venn_macro"]),
    "mal_mp3_to_mp1":         ("Mal_MP3→MP1",     COL["venn_mal"]),
    "neu_osm_priming_to_low": ("Neu_OSM_priming→low", COL["venn_neu"]),
}
N_MIN, FDR = 5, 0.05
TOP_NS = [50, 100, 200, 300]


def load_filtered(t):
    df = pd.read_csv(PERT / t / f"{t}_stats.csv", index_col=0)
    df = df[(df["Sig"] == 1) & (df["Shift_to_goal_end"] > 0)
            & (df["Goal_end_FDR"] < FDR) & (df["N_Detections"] >= N_MIN)].copy()
    return df.sort_values("Shift_to_goal_end", ascending=False).reset_index(drop=True)


filt = {t: load_filtered(t) for t in TRANSITIONS}
for t, d in filt.items():
    print(f"  {t}: {len(d)} candidates pass filter")

# also count *unfiltered* total candidates (Sig=1 only, any direction)
unfiltered_n = {}
for t in TRANSITIONS:
    df = pd.read_csv(PERT / t / f"{t}_stats.csv", index_col=0)
    unfiltered_n[t] = int((df["Sig"] == 1).sum())
print(f"unfiltered Sig=1 sizes: {unfiltered_n}")


def gset(d, n): return set(d.head(n)["Gene_name"])


def make_venn_with_bar(top_n, out_stem):
    a = gset(filt["macro_spp1_to_c1qc"], top_n)
    b = gset(filt["mal_mp3_to_mp1"], top_n)
    c = gset(filt["neu_osm_priming_to_low"], top_n)
    only_a, only_b, only_c = a - b - c, b - a - c, c - a - b
    ab, ac, bc = (a & b) - c, (a & c) - b, (b & c) - a
    abc = a & b & c

    fig = plt.figure(figsize=(5.6, 3.4))
    gs = GridSpec(1, 2, width_ratios=[3.5, 1.0], wspace=0.05, figure=fig)
    ax_v = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])

    subsets = (len(only_a), len(only_b), len(ab), len(only_c),
               len(ac), len(bc), len(abc))
    labels = [TRANSITIONS["macro_spp1_to_c1qc"][0],
              TRANSITIONS["mal_mp3_to_mp1"][0],
              TRANSITIONS["neu_osm_priming_to_low"][0]]
    cols_v = (TRANSITIONS["macro_spp1_to_c1qc"][1],
              TRANSITIONS["mal_mp3_to_mp1"][1],
              TRANSITIONS["neu_osm_priming_to_low"][1])
    v = venn3(subsets=subsets, set_labels=labels,
              set_colors=cols_v, alpha=0.55, ax=ax_v)
    venn3_circles(subsets=subsets, linewidth=0.6, ax=ax_v, color="black")
    for lbl in v.set_labels:
        if lbl is not None: lbl.set_fontsize(7)
    for lbl in v.subset_labels:
        if lbl is not None: lbl.set_fontsize(7)
    if abc:
        ax_v.text(0.5, -0.07, f"3-way: {', '.join(sorted(abc))}",
                  ha="center", fontsize=6, transform=ax_v.transAxes)
    add_subtitle(ax_v, f"top-{top_n} per transition · Sig=1 · FDR<{FDR}")

    # right barplot — per-transition pool size at this top-N
    bar_keys = list(TRANSITIONS.keys())
    bar_short = ["Macro", "Mal", "Neu"]
    bar_vals = [len(gset(filt[k], top_n)) for k in bar_keys]
    bar_cols = [TRANSITIONS[k][1] for k in bar_keys]
    bars = ax_b.barh(bar_short, bar_vals, color=bar_cols,
                     edgecolor="black", linewidth=0.4)
    for bar, val in zip(bars, bar_vals):
        ax_b.text(val + max(bar_vals) * 0.02, bar.get_y() + bar.get_height() / 2,
                  str(val), va="center", fontsize=6)
    style_axes(ax_b, xlabel="# candidate genes", ylabel="")
    ax_b.invert_yaxis()
    add_subtitle(ax_b, f"top-{top_n}")

    fig.tight_layout()
    save_panel(fig, out_stem)
    plt.close(fig)
    return abc, (ab, ac, bc)


pair_rows, shared_rows = [], []
for n in TOP_NS:
    abc, (ab, ac, bc) = make_venn_with_bar(n, OUT / f"8A_top{n}_venn")
    pair_rows.append({"top_N": n,
                      "macro_AND_mal": ";".join(sorted(ab)),
                      "macro_AND_neu": ";".join(sorted(ac)),
                      "mal_AND_neu":   ";".join(sorted(bc)),
                      "all_three":     ";".join(sorted(abc))})
    if abc:
        for g in sorted(abc):
            row = {"top_N": n, "Gene_name": g}
            for t, d in filt.items():
                hit = d[d["Gene_name"] == g]
                if len(hit):
                    row[f"{t}_rank"] = int(d.index[d["Gene_name"] == g][0]) + 1
                    row[f"{t}_shift"] = hit["Shift_to_goal_end"].iloc[0]
                    row[f"{t}_FDR"] = hit["Goal_end_FDR"].iloc[0]
            shared_rows.append(row)
pd.DataFrame(pair_rows).to_csv(OUT / "8A_pairwise.csv", index=False)
if shared_rows:
    pd.DataFrame(shared_rows).to_csv(OUT / "8A_shared_3way.csv", index=False)

# ≥2-transition pool at top-200
N_FOR_POOL = 200
a = gset(filt["macro_spp1_to_c1qc"], N_FOR_POOL)
b = gset(filt["mal_mp3_to_mp1"], N_FOR_POOL)
c = gset(filt["neu_osm_priming_to_low"], N_FOR_POOL)
pool = (a & b) | (a & c) | (b & c)
pool_rows = []
for g in sorted(pool):
    row = {"Gene_name": g}; hits = []
    for t, d in filt.items():
        sub = d[d["Gene_name"] == g]
        if len(sub):
            row[f"{t[:5]}_rank"]  = int(d.index[d["Gene_name"] == g][0]) + 1
            row[f"{t[:5]}_shift"] = sub["Shift_to_goal_end"].iloc[0]
            row[f"{t[:5]}_FDR"]   = sub["Goal_end_FDR"].iloc[0]
            hits.append(t)
    row["n_transitions_hit"] = len(hits)
    row["transitions_hit"] = ";".join(hits)
    for t, d in filt.items():
        sub = d[d["Gene_name"] == g]
        if len(sub):
            row["Ensembl_ID"] = sub["Ensembl_ID"].iloc[0]; break
    pool_rows.append(row)
pool_df = pd.DataFrame(pool_rows).sort_values(
    ["n_transitions_hit", "Gene_name"], ascending=[False, True])
pool_df.to_csv(OUT / "8A_candidate_pool.csv", index=False)
print(pool_df[["Gene_name", "Ensembl_ID", "n_transitions_hit",
               "transitions_hit"]].to_string(index=False))

# per-transition top table
rows = []
for t, d in filt.items():
    for n in TOP_NS:
        h = d.head(n).copy(); h["transition"] = t; h["top_N"] = n
        rows.append(h[["transition", "top_N", "Gene_name", "Ensembl_ID",
                       "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]])
pd.concat(rows, ignore_index=True).to_csv(OUT / "8A_per_transition_top.csv", index=False)

print("\nDONE.")
