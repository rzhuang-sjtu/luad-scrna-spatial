"""Step 9d: export per-section CSVs for R MISTy on Okamura cohort."""
from pathlib import Path
import pandas as pd, scanpy as sc

ROOT = Path("${DATA_ROOT}/ST/results/step09_okamura_validation")
OUT = ROOT / "misty_data"
OUT.mkdir(parents=True, exist_ok=True)

a = sc.read_h5ad(str(ROOT / "cohort_with_progeny.h5ad"))
print("loaded cohort:", a.shape)
abund = a.obsm["q05_cell_abundance"].copy()
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
print("[done]", OUT)
