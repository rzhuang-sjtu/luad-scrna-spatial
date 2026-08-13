"""Step 7c: clean MISTy heatmap — cell-type predictors × pathway targets, intra view."""
import os
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUT = Path("${DATA_ROOT}/ST/results/step07_misty")
agg = pd.read_csv(OUT / "aggregated_importance.csv")
print("rows:", len(agg), "views:", agg["view"].unique())

# keep only cell-type predictors (drop pathway-as-predictor rows)
ct_predictors = [p for p in agg["Predictor"].unique() if not p.startswith("progeny_")]
agg = agg[agg["Predictor"].isin(ct_predictors)].copy()
agg["target_name"] = agg["Target"].str.replace("progeny_", "").str.replace(".", "-", regex=False)

# Make heatmap per view
for view in agg["view"].unique():
    sub = agg[agg["view"] == view].pivot(index="Predictor", columns="target_name",
                                          values="mean_importance")
    if sub.empty: continue
    sub = sub.fillna(0)
    # order rows by mean magnitude
    sub = sub.reindex(sub.abs().mean(axis=1).sort_values(ascending=False).index)
    fig, ax = plt.subplots(figsize=(max(7, 0.8*sub.shape[1]+3), max(8, 0.4*sub.shape[0]+1)))
    sns.heatmap(sub, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                cbar_kws={"label": "mean MISTy importance"}, ax=ax)
    ax.set_title(f"MISTy importance ({view}): cell type → PROGENy pathway")
    ax.set_xlabel("PROGENy pathway")
    ax.set_ylabel("cell type")
    fig.tight_layout()
    out_png = OUT / f"clean_heatmap_{view.replace('.','_')}.png"
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved", out_png)
print("[done]")
