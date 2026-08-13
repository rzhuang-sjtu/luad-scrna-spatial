"""Fig 8H/I/J/K/L + S11A-C — ST candidate-gene panels.

Panel layout per checklist:
  8H : candidate-gene co-expression UMAP on ST spots, colored by 3 candidate genes (3 sub-panels).
  8I : 1 high-MP4 (R-surrogate) section × 3 candidate genes (1×3 spatial grid).
  8J : 2nd  high-MP4 section × 3 genes.
  8K : 1 high-MP3 (NR-surrogate) section × 3 genes.
  8L : 2nd high-MP3 section × 3 genes.
  S11A-C : 3 additional sections (mixed/extra) × 3 genes.

R/NR surrogate definition (LUAD): per-section mean(MP4_score) vs mean(MP3_score).
  Top-2 by MP4 → R-surrogate sections (8I, 8J).
  Top-2 by MP3 → NR-surrogate sections (8K, 8L).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from scipy import stats
import scipy.sparse as sp

import sys
from pathlib import Path
# fig8_style lives with the figure it styles; make it importable when this
# script is run from its own directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "plotting" / "fig8"))
from fig8_style import setup, COL, CMAP_ST, save_panel, sig_stars, add_subtitle, style_axes
setup()

ROI_H5 = Path("${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad")
HE_ROOT = Path("${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530")
# tissue_hires_scalef in scalefactors.json (uniform across sections)
HIRES_SCALEF = 0.11586143
CAND_FILE = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn/8A_candidate_pool.csv")
OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8H_st")
OUT.mkdir(parents=True, exist_ok=True)

cand = pd.read_csv(CAND_FILE)
ALL_CANDS = list(cand["Gene_name"])
TOP3 = ["HSP90AB1", "SRSF9", "NDUFB2"]

print("loading h5ad ...")
ad = sc.read_h5ad(ROI_H5)
sample_col = "sample"
print(f"  ad: {ad.shape}; samples: {ad.obs[sample_col].nunique()}")

# section-level MP scores → R/NR surrogate
agg = ad.obs.groupby(sample_col)[["MP3_score", "MP4_score"]].mean()
print("\nper-section mean MP scores:")
print(agg.to_string())

# R surrogate = top-2 by MP4 score; NR surrogate = top-2 by MP3 score
sec_R  = agg.sort_values("MP4_score", ascending=False).head(2).index.tolist()
sec_NR = agg.sort_values("MP3_score", ascending=False).head(2).index.tolist()
# extra (S11A-C): three sections not in R/NR
used = set(sec_R) | set(sec_NR)
sec_extra = [s for s in agg.index if s not in used][:3]
print(f"\nR surrogate (high MP4, 8I/J): {sec_R}")
print(f"NR surrogate (high MP3, 8K/L): {sec_NR}")
print(f"S11A-C extra: {sec_extra}")

def load_he(section):
    """Return (img, scalef) where img is the tissue_hires_image and scalef is the
    factor that converts full-res spatial coords to hires-image pixel coords."""
    p = HE_ROOT / f"{section}-spatial" / "tissue_hires_image.png"
    if not p.exists():
        return None, None
    return mpimg.imread(str(p)), HIRES_SCALEF


def plot_section_3genes(section, genes, out_stem, subtitle_prefix=""):
    mask = ad.obs[sample_col].astype(str) == section
    sub = ad[mask]
    if "spatial" in sub.obsm:
        xy_full = sub.obsm["spatial"]
    elif {"x", "y"}.issubset(sub.obs.columns):
        xy_full = sub.obs[["x", "y"]].values
    else:
        print(f"  WARN: no spatial coord for {section}")
        return

    he, scalef = load_he(section)
    if he is None:
        print(f"  WARN: no H&E for {section}")
        scalef = 1.0
    # spatial coords in cohort_with_roi.h5ad are stored in full-res pixel space;
    # multiply by hires scalef to align with tissue_hires_image.png
    xy = xy_full * scalef

    fig, axes = plt.subplots(1, len(genes), figsize=(2.8 * len(genes), 2.9),
                             squeeze=False)
    for j, g in enumerate(genes):
        ax = axes[0, j]
        # 1) H&E underlay
        if he is not None:
            ax.imshow(he, origin="upper", interpolation="bilinear", alpha=0.85)
        if g not in sub.var_names:
            ax.text(0.5, 0.5, f"{g}\nnot in ST", ha="center",
                    transform=ax.transAxes, fontsize=7)
            ax.set_xticks([]); ax.set_yticks([]); continue
        # 2) spot scatter colored by gene — blue→red sequential
        e = sub[:, g].X
        if hasattr(e, "toarray"): e = e.toarray().flatten()
        else: e = np.asarray(e).flatten()
        h = ax.scatter(xy[:, 0], xy[:, 1], c=e, s=3.5, cmap=CMAP_ST,
                       edgecolor="none", alpha=0.78, rasterized=True)
        if he is not None:
            xmin, xmax = xy[:, 0].min(), xy[:, 0].max()
            ymin, ymax = xy[:, 1].min(), xy[:, 1].max()
            mx = (xmax - xmin) * 0.04; my = (ymax - ymin) * 0.04
            ax.set_xlim(max(0, xmin - mx), min(he.shape[1], xmax + mx))
            ax.set_ylim(min(he.shape[0], ymax + my), max(0, ymin - my))
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        for sp_ in ax.spines.values(): sp_.set_visible(False)
        # gene name above each subplot
        add_subtitle(ax, g)
        cb = plt.colorbar(h, ax=ax, fraction=0.035, pad=0.015, shrink=0.85)
        cb.ax.tick_params(labelsize=6)
        cb.outline.set_linewidth(0.4)
    # section info as a footer below the row of subplots
    fig.tight_layout()
    if subtitle_prefix:
        fig.subplots_adjust(bottom=0.10)
        fig.text(0.5, 0.02, subtitle_prefix, fontsize=7, style="italic",
                 ha="center", va="bottom", transform=fig.transFigure)
    save_panel(fig, out_stem)
    plt.close(fig)


print("\n=== 8I/J — R surrogate (high MP4) ===")
for letter, sec in zip(["I", "J"], sec_R):
    mp3 = agg.loc[sec, "MP3_score"]; mp4 = agg.loc[sec, "MP4_score"]
    plot_section_3genes(sec, TOP3, OUT / f"8{letter}_spatial_R_{sec}",
                        subtitle_prefix=f"{sec} (R surrogate · MP4={mp4:.2f}, MP3={mp3:.2f})")
    print(f"  8{letter} → {sec}")

print("\n=== 8K/L — NR surrogate (high MP3) ===")
for letter, sec in zip(["K", "L"], sec_NR):
    mp3 = agg.loc[sec, "MP3_score"]; mp4 = agg.loc[sec, "MP4_score"]
    plot_section_3genes(sec, TOP3, OUT / f"8{letter}_spatial_NR_{sec}",
                        subtitle_prefix=f"{sec} (NR surrogate · MP3={mp3:.2f}, MP4={mp4:.2f})")
    print(f"  8{letter} → {sec}")

print("\n=== S11A-C — extra sections ===")
for letter, sec in zip(["A", "B", "C"], sec_extra):
    mp3 = agg.loc[sec, "MP3_score"]; mp4 = agg.loc[sec, "MP4_score"]
    plot_section_3genes(sec, TOP3, OUT / f"S11{letter}_spatial_{sec}",
                        subtitle_prefix=f"{sec} (MP3={mp3:.2f}, MP4={mp4:.2f})")
    print(f"  S11{letter} → {sec}")

# Approach: subsample 8000 spots, normalize+log+scale TOP3 only, run UMAP on this 3D
# space (n_components=2), color each sub-panel by one gene.
print("\n=== 8H — candidate-gene co-expression UMAP ===")

import scanpy as sc_mod
rng = np.random.default_rng(0)
n_max = 8000
ad_idx = np.arange(ad.n_obs)
if ad.n_obs > n_max:
    keep_idx = rng.choice(ad_idx, size=n_max, replace=False)
else:
    keep_idx = ad_idx
sub = ad[keep_idx, :].copy()
# get expr matrix for the candidate set
cand_in_st = [g for g in ALL_CANDS if g in sub.var_names]
expr_X = sub[:, cand_in_st].X
if hasattr(expr_X, "toarray"): expr_X = expr_X.toarray()
expr_df = pd.DataFrame(expr_X, columns=cand_in_st, index=sub.obs_names)
# scale per-gene
expr_z = (expr_df - expr_df.mean()) / expr_df.std().replace(0, 1)

# UMAP on candidate-gene space
import umap
um = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0).fit_transform(
    expr_z.fillna(0).values)
print(f"  UMAP done: {um.shape}")

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4), squeeze=False)
for i, g in enumerate(TOP3):
    ax = axes[0, i]
    if g not in cand_in_st:
        ax.axis("off"); continue
    h = ax.scatter(um[:, 0], um[:, 1], c=expr_df[g].values, s=1.2,
                   cmap="viridis", edgecolor="none", alpha=0.9, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_visible(False)
    add_subtitle(ax, g)
    ax.set_xlabel("UMAP-1", fontsize=7); ax.set_ylabel("UMAP-2" if i == 0 else "", fontsize=7)
    cb = plt.colorbar(h, ax=ax, fraction=0.04, pad=0.02)
    cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.4)
fig.text(0.01, 0.98, f"Candidate-gene co-expression UMAP · n={len(keep_idx)} ST spots",
         fontsize=7, style="italic", va="top")
fig.tight_layout()
save_panel(fig, OUT / "8H_candidate_umap")
plt.close(fig)
print("  wrote 8H_candidate_umap")

roi_col = "roi" if "roi" in ad.obs.columns else next(c for c in ad.obs.columns if "roi" in c.lower())
expr_full = ad[:, [g for g in TOP3 if g in ad.var_names]].X
if hasattr(expr_full, "toarray"): expr_full = expr_full.toarray()
ed = pd.DataFrame(expr_full, columns=[g for g in TOP3 if g in ad.var_names])
ed["roi"] = ad.obs[roi_col].astype(bool).values
stat_rows = []
for g in TOP3:
    if g not in ed.columns: continue
    a = ed.loc[ed["roi"], g]; b = ed.loc[~ed["roi"], g]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    stat_rows.append({"gene": g, "roi_mean": a.mean(), "nonroi_mean": b.mean(),
                      "delta": a.mean() - b.mean(), "wilcoxon_p": p})
pd.DataFrame(stat_rows).to_csv(OUT / "8I_roi_vs_nonroi_stats.csv", index=False)

# section table
agg.to_csv(OUT / "8H_section_MP_scores.csv")
print("\nDONE.")
