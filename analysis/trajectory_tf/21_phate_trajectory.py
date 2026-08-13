"""Step 21 (Fig 3B): PHATE 2D embedding for Monocle-like Component_1/Component_2 axes.

PHATE produces smoother, more trajectory-like 2D layouts than UMAP — closer to
Monocle3 DDRTree visual. Run on X_pca_mal_harmony (50 PCs) for 45k malignant
cells, embed to 2D, export.

Outputs (overwriting old monocle_trajectory.csv.gz):
  ${WORK_ROOT}/luad_figures/fig2/monocle_trajectory.csv.gz
    barcode, dataset, Component_1, Component_2, pseudotime, pseudotime_rank, dominant_MP
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
PT_CSV = Path.home()/"luad/results/step10b_pseudotime.csv.gz"
FIG2 = Path("${WORK_ROOT}/luad_figures/fig2")


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    pca = a.obsm["X_pca_mal_harmony"]
    log(f"  cells={a.n_obs}, harmony PCs={pca.shape}")

    log("running PHATE on Harmony PCs")
    import phate
    op = phate.PHATE(n_components=2, n_jobs=8, random_state=0,
                     knn=15, decay=40, t="auto", verbose=1)
    embed = op.fit_transform(pca)
    log(f"  PHATE embed: {embed.shape}, range "
        f"X=[{embed[:,0].min():.3f},{embed[:,0].max():.3f}] "
        f"Y=[{embed[:,1].min():.3f},{embed[:,1].max():.3f}]")

    # Pseudotime + rank-transform
    pt = pd.read_csv(PT_CSV)[["barcode","pseudotime"]]
    pt_lookup = dict(zip(pt["barcode"], pt["pseudotime"]))
    pseudotime = np.array([pt_lookup.get(bc, np.nan) for bc in a.obs.index])
    valid = np.isfinite(pseudotime)
    pt_rank = np.full_like(pseudotime, np.nan, dtype=np.float64)
    ranks = pd.Series(pseudotime[valid]).rank(pct=True).values
    pt_rank[valid] = ranks
    log(f"  pseudotime: orig range [{pseudotime[valid].min():.4f},"
        f"{pseudotime[valid].max():.4f}]; rank-transformed to [0,1]")

    out = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "Component_1": embed[:, 0],
        "Component_2": embed[:, 1],
        "pseudotime": pseudotime,
        "pseudotime_rank": pt_rank,
        "dominant_MP": a.obs["dominant_MP"].astype(str).values,
    })
    out.to_csv(FIG2/"monocle_trajectory.csv.gz", index=False, compression="gzip")
    log(f"  monocle_trajectory.csv.gz overwritten ({out.shape})")
    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
