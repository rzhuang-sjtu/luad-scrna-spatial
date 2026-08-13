"""Step 10b (Fig 2G-H): Diffusion pseudotime on LUAD malignant cells.

Pipeline:
  1. Load luad_malignant_scored.h5ad (45k cells, has X_pca_mal_harmony + X_umap_mal + MP scores).
  2. Rebuild neighbors on X_pca_mal_harmony (stored neighbors were on joint Harmony).
  3. sc.tl.diffmap → sc.tl.dpt(n_branchings=0, 1).
  4. Root cell = cell with max MP4 score (AT2-like = most differentiated).
  5. Export UMAP + pseudotime + MP scores; bin-wise MP dominance density.

Outputs:
  - ~/luad/results/step10b_pseudotime.csv.gz
  - ${WORK_ROOT}/luad_figures/fig2/pseudotime_umap.csv.gz
  - ${WORK_ROOT}/luad_figures/fig2/pseudotime_mp_density.csv
  - ${WORK_ROOT}/luad_figures/fig2/pseudotime_summary.md
"""

from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home() / "luad/data/processed/luad_malignant_scored.h5ad"
OUTCSV = Path.home() / "luad/results/step10b_pseudotime.csv.gz"
FIGDIR = Path("${WORK_ROOT}/luad_figures/fig2")
FIGDIR.mkdir(parents=True, exist_ok=True)
VALID_MPS = ["MP1", "MP2", "MP3", "MP4"]
N_BINS = 30


def main() -> None:
    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] loading {IN}")
    a = sc.read_h5ad(IN)
    print(f"  shape={a.shape}  PCs={a.obsm['X_pca_mal_harmony'].shape}")

    # Rebuild neighbors on malignant Harmony rep (stored uses joint harmony)
    print(f"[{time.strftime('%H:%M:%S')}] rebuilding neighbors on X_pca_mal_harmony ...")
    sc.pp.neighbors(a, use_rep="X_pca_mal_harmony", n_neighbors=30, random_state=0)
    print(f"  neighbors built")

    print(f"[{time.strftime('%H:%M:%S')}] diffusion map (15 comp)...")
    sc.tl.diffmap(a, n_comps=15, random_state=0)

    # Root cell: MP4 (AT2-like) most differentiated cell
    mp4 = a.obs["MP4_score"].values
    iroot = int(np.argmax(mp4))
    root_bc = a.obs.index[iroot]
    print(f"  iroot={iroot} ({root_bc}), MP4_score={mp4[iroot]:.3f}, "
          f"dominant_MP={a.obs['dominant_MP'].iloc[iroot]}")
    a.uns["iroot"] = iroot

    print(f"[{time.strftime('%H:%M:%S')}] DPT (n_branchings=1)...")
    try:
        sc.tl.dpt(a, n_branchings=1)
    except Exception as e:
        print(f"  DPT with 1 branching failed ({e}); retrying with n_branchings=0")
        sc.tl.dpt(a, n_branchings=0)

    # Collect result
    pt = a.obs["dpt_pseudotime"].values
    print(f"  pseudotime range: [{np.nanmin(pt):.3f}, {np.nanmax(pt):.3f}]")
    print(f"  infinite/nan cells: {int((~np.isfinite(pt)).sum())} -> masking")

    umap = a.obsm.get("X_umap_mal", a.obsm.get("X_umap"))
    umap_key = "X_umap_mal" if "X_umap_mal" in a.obsm else "X_umap"

    df = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "patient_id": a.obs["patient_id"].astype(str).values,
        "tissue_type": a.obs["tissue_type"].astype(str).values,
        "UMAP1": umap[:, 0],
        "UMAP2": umap[:, 1],
        "pseudotime": pt,
        "dpt_groups": a.obs["dpt_groups"].astype(str).values if "dpt_groups" in a.obs else "-",
        "dominant_MP": a.obs["dominant_MP"].astype(str).values,
        "MP1_score": a.obs["MP1_score"].values,
        "MP2_score": a.obs["MP2_score"].values,
        "MP3_score": a.obs["MP3_score"].values,
        "MP4_score": a.obs["MP4_score"].values,
    })
    df.to_csv(OUTCSV, index=False, compression="gzip")

    # Fig 2G: UMAP + pseudotime + MP scores
    df[["barcode", "UMAP1", "UMAP2", "pseudotime", "dominant_MP",
        "MP1_score", "MP2_score", "MP3_score", "MP4_score"]].to_csv(
        FIGDIR / "pseudotime_umap.csv.gz", index=False, compression="gzip"
    )

    # Fig 2H: bin pseudotime → dominant MP proportion
    df_valid = df[np.isfinite(df["pseudotime"])].copy()
    df_valid["pt_bin"] = pd.qcut(df_valid["pseudotime"], q=N_BINS, labels=False,
                                  duplicates="drop")
    ct = pd.crosstab(df_valid["pt_bin"], df_valid["dominant_MP"])
    ct = ct.reindex(columns=VALID_MPS + (["MP5"] if "MP5" in ct.columns else []),
                    fill_value=0)
    # Include mean pseudotime per bin for plotting on x-axis
    bin_mean = df_valid.groupby("pt_bin")["pseudotime"].mean()
    density = ct.div(ct.sum(axis=1), axis=0).round(5)
    density.insert(0, "n_cells", ct.sum(axis=1))
    density.insert(0, "mean_pseudotime", bin_mean)
    density.to_csv(FIGDIR / "pseudotime_mp_density.csv")

    # Summary md
    print(f"[{time.strftime('%H:%M:%S')}] writing summary ...")
    with open(FIGDIR / "pseudotime_summary.md", "w") as f:
        f.write("# LUAD Malignant — Diffusion Pseudotime (Fig 2G-H)\n\n")
        f.write(f"- N cells: {len(df)}\n")
        f.write(f"- Valid pseudotime cells: {len(df_valid)}\n")
        f.write(f"- Root cell: {root_bc} (MP4 dominant, MP4_score={mp4[iroot]:.3f})\n")
        f.write(f"- dpt_groups count: "
                f"{a.obs['dpt_groups'].value_counts().to_dict() if 'dpt_groups' in a.obs else 'none'}\n\n")
        f.write("## Pseudotime distribution by dominant MP\n\n")
        stats = (df_valid.groupby("dominant_MP")["pseudotime"]
                 .agg(["mean", "median", "std", "count"]).round(3))
        f.write(stats.to_markdown() + "\n\n")
        f.write("## Pseudotime bin × dominant-MP proportion (first + last 5 bins)\n\n")
        f.write(density.head(5).round(3).to_markdown() + "\n\n...\n")
        f.write(density.tail(5).round(3).to_markdown() + "\n")

    # Save enriched h5ad (minimal additions to obs + obsm diffmap)
    # Keep in separate file to avoid huge rewrite
    add = a.obs[["dpt_pseudotime"]].copy()
    if "dpt_groups" in a.obs:
        add["dpt_groups"] = a.obs["dpt_groups"]
    add.to_csv(Path.home() / "luad/results/step10b_pseudotime_obs.csv.gz",
               compression="gzip")

    print(f"\n=== Pseudotime x dominant MP ===")
    print(stats)
    print(f"\n[{time.strftime('%H:%M:%S')}] done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
