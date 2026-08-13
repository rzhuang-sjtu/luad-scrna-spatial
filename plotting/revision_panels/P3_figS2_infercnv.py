#!/usr/bin/env python
"""New Fig. S2 panel: inferCNV cross-validation of the CopyKAT malignant call.

Revision analysis:. The submitted Limitations conceded that malignant-cell
identification was never checked against an independent CNV method; it now has
been, and the manuscript states the result with numbers, so the figure has to
carry it.

Drawn with the manuscript's own fig8_style so it sits beside the published
Fig. S2 panels without restyling.

Output: revision_panels/FigS2_new__infercnv_five_tiers.{pdf,png}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "fig8"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fig8_style import setup, COL, save_panel, add_subtitle, style_axes

setup()
plt.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.default": "regular",
})

RES = Path("${PROJECT_ROOT}/results/infercnv")
OUT = Path("${WORK_ROOT}/revision_panels")

d = pd.read_csv(RES / "cell_level.csv.gz")
auroc = pd.read_csv(RES / "auroc.csv").set_index("metric").loc["frac_altered"]

# Five tiers, ordered from the diploid reference to the malignant call.
tier = np.where(
    d.group.eq("reference"), "Reference cells\n(T/NK, B, myeloid)",
    np.where(d.celltype_coarse.ne("Epithelial"),
             "Non-epithelial cells\nin the observation set",
             np.where(d.malignant.eq("Non-malignant"),
                      "Epithelial,\nCopyKAT non-malignant",
                      np.where(d.malignant.eq("Uncertain"),
                               "Epithelial,\nCopyKAT uncertain",
                               "Epithelial,\nCopyKAT malignant"))))
d = d.assign(tier=tier)
order = ["Reference cells\n(T/NK, B, myeloid)",
         "Non-epithelial cells\nin the observation set",
         "Epithelial,\nCopyKAT non-malignant",
         "Epithelial,\nCopyKAT uncertain",
         "Epithelial,\nCopyKAT malignant"]

fig, ax = plt.subplots(figsize=(5.0, 2.6))
data = [d.loc[d.tier.eq(t), "frac_altered"].values for t in order]
bp = ax.boxplot(data, vert=False, widths=0.62, showfliers=False,
                patch_artist=True, medianprops=dict(color="black", lw=0.9),
                whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                boxprops=dict(lw=0.5))
# One accent: only the malignant tier is coloured, everything else neutral.
for i, box in enumerate(bp["boxes"]):
    box.set_facecolor(COL["tumor"] if i == len(order) - 1 else "#D9D9D9")
    box.set_edgecolor("black")

for i, t in enumerate(order, start=1):
    v = d.loc[d.tier.eq(t), "frac_altered"]
    ax.text(1.005, i, f"n = {len(v):,}", transform=ax.get_yaxis_transform(),
            va="center", ha="left", fontsize=6, color="#4D4D4D")

ax.set_yticks(range(1, len(order) + 1))
ax.set_yticklabels(order, fontsize=6.5)
ax.set_xlim(0, 1.0)
style_axes(ax, xlabel="Fraction of the genome altered (inferCNV)", ylabel="")
add_subtitle(ax, "161,765 cells from 89 patients with sufficient diploid "
                 "reference · boxes, median and IQR")
ax.text(0.985, 0.06,
        "malignant vs non-malignant epithelium\n"
        "AUROC = %.3f (95%% CI %.3f–%.3f)\n77,796 cells, 74 patients"
        % (auroc.auroc_overall, auroc.ci_lo, auroc.ci_hi),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6,
        color="#4D4D4D")
fig.tight_layout()
save_panel(fig, str(OUT / "FigS2_new__infercnv_five_tiers"))
plt.close(fig)

med = d.groupby("tier")["frac_altered"].median().reindex(order)
print("median fraction altered per tier:")
for t, v in med.items():
    print(f"  {t.replace(chr(10), ' '):48s} {v:.4f}")
print(f"AUROC {auroc.auroc_overall:.3f} "
      f"({auroc.ci_lo:.3f}-{auroc.ci_hi:.3f}), "
      f"{int(auroc.n_cells):,} cells, {int(auroc.n_patients_evaluable)} patients")
print("written", OUT / "FigS2_new__infercnv_five_tiers")
