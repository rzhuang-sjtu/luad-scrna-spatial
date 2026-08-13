"""Step 6b: re-plot COMMOT results from saved section h5ads, plus summary."""
import os, time
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT    = Path("${DATA_ROOT}/ST/results/step06_commot")
SEC    = OUT / "section_h5ad"
PLOTS  = OUT / "spatial_plots"
PLOTS.mkdir(exist_ok=True, parents=True)

PATHWAYS = ["OSM", "IL1"]
LOG = OUT / "run.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")


samples = sorted([p.stem for p in SEC.glob("*.h5ad")])
log(f"replotting {len(samples)} sections")

rows = []
for s in samples:
    a = sc.read_h5ad(str(SEC / f"{s}.h5ad"))
    a.var_names_make_unique()
    sender_df = a.obsm["commot-cellchat-sum-sender"]
    recv_df   = a.obsm["commot-cellchat-sum-receiver"]
    # Sum sender / receiver across all LR pairs of a pathway
    for p in PATHWAYS:
        s_cols = [c for c in sender_df.columns if c.startswith(f"s-") and f"-{p}-" in f"-{c}-".replace("s-","")[:200] or
                                                  c == f"s-{p}"]
        r_cols = [c for c in recv_df.columns   if c.startswith(f"r-") and f"-{p}-" in f"-{c}-".replace("r-","")[:200] or
                                                  c == f"r-{p}"]
        # safer match: pathway either appears as suffix-then-LR or exact column f"s-{p}"
        s_cols_strict = [c for c in sender_df.columns if c == f"s-{p}" or c.startswith(f"s-{p}-") or
                          (c.startswith("s-") and c.split("-",2)[1] == p)]
        r_cols_strict = [c for c in recv_df.columns   if c == f"r-{p}" or c.startswith(f"r-{p}-") or
                          (c.startswith("r-") and c.split("-",2)[1] == p)]
        # use strict
        s_cols, r_cols = s_cols_strict, r_cols_strict
        a.obs[f"s_{p}"] = sender_df[s_cols].sum(axis=1).values if s_cols else 0.0
        a.obs[f"r_{p}"] = recv_df[r_cols].sum(axis=1).values    if r_cols else 0.0
        a.obs[f"total_{p}"] = a.obs[f"s_{p}"].fillna(0) + a.obs[f"r_{p}"].fillna(0)

    rows.append({
        "sample": s,
        "n_spots": a.n_obs,
        "OSM_send_mean": float(a.obs["s_OSM"].mean()),
        "OSM_recv_mean": float(a.obs["r_OSM"].mean()),
        "OSM_send_pct_active": float((a.obs["s_OSM"] > 0).mean()),
        "OSM_recv_pct_active": float((a.obs["r_OSM"] > 0).mean()),
        "IL1_send_mean": float(a.obs["s_IL1"].mean()),
        "IL1_recv_mean": float(a.obs["r_IL1"].mean()),
        "IL1_send_pct_active": float((a.obs["s_IL1"] > 0).mean()),
        "IL1_recv_pct_active": float((a.obs["r_IL1"] > 0).mean()),
    })

    # Plot 2x3 grid: per pathway sender / receiver / total
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for i, p in enumerate(PATHWAYS):
        for j, kind in enumerate(["s", "r", "total"]):
            col = f"{kind}_{p}" if kind != "total" else f"total_{p}"
            title = f"{p} {'sender' if kind=='s' else 'receiver' if kind=='r' else 'total'}"
            ax = axes[i, j]
            try:
                sc.pl.spatial(a, color=col, library_id=s, ax=ax, show=False,
                              cmap="magma", size=1.4, frameon=False, title=title, colorbar_loc="right")
            except Exception as e:
                ax.set_title(f"{title} (err)")
                ax.text(0.5, 0.5, f"{type(e).__name__}", ha="center", va="center")
    fig.suptitle(f"{s}  COMMOT (CellChat)  OSM + IL1", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / f"{s}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # Save updated h5ad with extracted obs columns
    a.write_h5ad(str(SEC / f"{s}.h5ad"), compression="gzip")
    log(f"  {s}: OSM s_mean={rows[-1]['OSM_send_mean']:.3f} r_mean={rows[-1]['OSM_recv_mean']:.3f}; "
        f"IL1 s_mean={rows[-1]['IL1_send_mean']:.3f} r_mean={rows[-1]['IL1_recv_mean']:.3f}")

df = pd.DataFrame(rows)
df.to_csv(OUT / "per_sample_pathway_summary.csv", index=False)
log("[done]")
print(df.round(3).to_string(index=False))
