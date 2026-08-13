"""Step 9e: clean MISTy heatmap for Okamura cohort."""
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from pathlib import Path
OUT = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/misty")
agg = pd.read_csv(OUT / "aggregated_importance.csv")
agg = agg[~agg["Predictor"].str.startswith("progeny_")].copy()
agg["target_name"] = agg["Target"].str.replace("progeny_", "").str.replace(".", "-", regex=False)
for view in agg["view"].unique():
    sub = agg[agg["view"] == view].pivot(index="Predictor", columns="target_name",
                                          values="mean_importance").fillna(0)
    sub = sub.reindex(sub.abs().mean(axis=1).sort_values(ascending=False).index)
    fig, ax = plt.subplots(figsize=(max(7, 0.8*sub.shape[1]+3), max(8, 0.4*sub.shape[0]+1)))
    sns.heatmap(sub, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                cbar_kws={"label": "mean MISTy importance"}, ax=ax)
    ax.set_title(f"Okamura MISTy ({view}): cell type -> PROGENy pathway")
    ax.set_xlabel("PROGENy pathway"); ax.set_ylabel("cell type")
    fig.tight_layout()
    out_png = OUT / f"clean_heatmap_{view.replace('.', '_')}.png"
    fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("saved", out_png)
print("[done]")
