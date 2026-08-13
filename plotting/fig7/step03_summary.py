"""Quick summary of per-sample mean q05 cell abundances after Step 3."""
import pandas as pd
from pathlib import Path

OUT = Path("${DATA_ROOT}/ST/results/step03_deconvolution")
d = pd.read_csv(OUT / "per_sample_mean_q05_abundance.csv", index_col=0)
d.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in d.columns]

groups = {
    "Neu subtypes": [c for c in d.columns if c.startswith("Neu_")],
    "Macro subtypes": [c for c in d.columns if c.startswith("Macro_")],
    "DC / Mono": [c for c in d.columns if c.startswith("cDC") or c.startswith("Mono") or c == "pDC"],
    "Malignant + stromal + lymphoid": ["Malignant", "Epithelial", "Fibroblast", "Endothelial", "T_NK", "B", "Plasma", "Mast"],
}
for title, cols in groups.items():
    cols = [c for c in cols if c in d.columns]
    print(f"\n=== {title} ===")
    print(d[cols].round(3).to_string())

# Top cell type per sample
print("\n=== top-3 cell types per sample (by mean q05 abundance) ===")
for s in d.index:
    row = d.loc[s].sort_values(ascending=False)
    top3 = ", ".join(f"{c}={v:.2f}" for c, v in row.head(3).items())
    print(f"   {s}: {top3}")
