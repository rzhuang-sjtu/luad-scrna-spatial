"""
Dump per-spot SEC61G/SRSF9/ANGPTL4 expression for both ST cohorts so
panels_R_part2.R can re-render the atlases with the corrected (Takano)
labels and shared-legend design.

Inputs:
  E-MTAB-13530: ${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad
  Takano (Okamura key) : ${DATA_ROOT}/ST/results/step09_okamura_validation/cohort.h5ad
                          (the original input cohort — c2l outputs were lost
                           in the May-7 crash; raw expression is sufficient
                           for the atlas because we only need per-spot expr)

Output:
  ${WORK_ROOT}/luad_figures/fig8/v2_500/data/8I_spot_long.csv
    Schema: cohort, sample, panel, kind, gene, spatial1, spatial2, expr
  Backup of previous CSV at 8I_spot_long.csv.bak
"""
from __future__ import annotations
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import scanpy as sc

GENES = ["SEC61G", "SRSF9", "ANGPTL4"]

E_PATH = Path("${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad")
O_PATH = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/cohort.h5ad")
PANEL_CSV = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data/8I_panel_assignments.csv")
OUT = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data/8I_spot_long.csv")

# Restrict Takano to the 8 fresh-frozen sections used in figs (the
# 8 FFPE were prepared but never deconvolved).
TAKANO_FF = {f"LUAD_No_{n}" for n in (1, 2, 3, 4, 5, 14, 16, 17)}

panel_assignments = pd.read_csv(PANEL_CSV)
print(f"panel_assignments rows: {len(panel_assignments)}")
panel_lookup = panel_assignments.set_index(["sample"]).to_dict(orient="index")

def dump(path, cohort_label, sample_filter=None):
    print(f"\n=== {cohort_label}: {path} ===")
    ad = sc.read_h5ad(path)
    if sample_filter is not None:
        keep = ad.obs["sample"].isin(sample_filter)
        ad = ad[keep].copy()
    print(f"  shape: {ad.shape}; samples: {sorted(ad.obs['sample'].unique())}")
    rows = []
    spatial = ad.obsm["spatial"]    # (n, 2)
    samples = ad.obs["sample"].values
    for g in GENES:
        if g not in ad.var_names:
            print(f"  [SKIP] {g} not in var_names"); continue
        e = ad[:, g].X
        if hasattr(e, "toarray"):
            e = e.toarray().flatten()
        else:
            e = np.asarray(e).flatten()
        for i in range(ad.n_obs):
            s = samples[i]
            info = panel_lookup.get(s, {})
            rows.append({
                "cohort":   cohort_label,
                "sample":   s,
                "panel":    info.get("panel", ""),
                "kind":     info.get("kind", ""),
                "gene":     g,
                "spatial1": float(spatial[i, 0]),
                "spatial2": float(spatial[i, 1]),
                "expr":     float(e[i]),
            })
    return rows

all_rows = []
all_rows += dump(E_PATH, "E-MTAB-13530")
all_rows += dump(O_PATH, "Okamura", sample_filter=TAKANO_FF)
df = pd.DataFrame(all_rows)
print(f"\ntotal rows: {len(df)}")

# Backup + write
if OUT.exists():
    shutil.copy(OUT, OUT.with_suffix(".csv.bak"))
df.to_csv(OUT, index=False)
print(f"wrote {OUT}")
print(df.groupby(["cohort","sample"]).size().rename("n_spots").to_string())
