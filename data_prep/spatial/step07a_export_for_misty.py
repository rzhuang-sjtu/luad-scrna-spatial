"""Step 7a: export per-section c2l abundance + PROGENy scores + spatial coords for R MISTy."""
import os, time
from pathlib import Path
import pandas as pd, numpy as np, scanpy as sc

COHORT = Path("${DATA_ROOT}/ST/results/step05_progeny/cohort_with_progeny.h5ad")
OUT = Path("${DATA_ROOT}/ST/results/step07_misty/data")
OUT.mkdir(parents=True, exist_ok=True)

a = sc.read_h5ad(str(COHORT))
print("loaded cohort:", a.shape)

abund = a.obsm["q05_cell_abundance"].copy()                     # cell-type abundances (intra view target features)
progeny_cols = [c for c in a.obs.columns if c.startswith("progeny_")]
mp_cols = [c for c in a.obs.columns if c.endswith("_score") and c.startswith("MP")]
print("intra cell types:", abund.shape, "| progeny cols:", len(progeny_cols), "| MP cols:", len(mp_cols))

samples = sorted(a.obs["sample"].unique().tolist())
for s in samples:
    mask = a.obs["sample"].values == s
    sub_idx = a.obs_names[mask]
    coords = a.obsm["spatial"][mask]
    abund_s   = abund.loc[sub_idx]
    progeny_s = a.obs.loc[sub_idx, progeny_cols]
    mp_s      = a.obs.loc[sub_idx, mp_cols]
    coords_df = pd.DataFrame(coords, index=sub_idx, columns=["x","y"])
    abund_s.to_csv(OUT / f"{s}_intra.csv")
    progeny_s.to_csv(OUT / f"{s}_progeny.csv")
    mp_s.to_csv(OUT / f"{s}_mp.csv")
    coords_df.to_csv(OUT / f"{s}_coords.csv")
    print(f"  {s}: spots={len(sub_idx)}")
print("[done] exported", len(samples), "sections to", OUT)
