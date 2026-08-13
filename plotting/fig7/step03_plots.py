"""
Step 3.3: Quick spatial sanity plots for 12 LUAD tumor sections.

Each section -> a 2x4 multi-panel figure with cell2location q05 abundance
heatmaps for: Malignant, Macro_SPP1, Neu_Inflammatory, Neu_OSM_priming,
Neu_OSM_low, Fibroblast, Endothelial.

Outputs: ${DATA_ROOT}/ST/results/step03_deconvolution/spatial_plots/<sample>.png
"""
from __future__ import annotations
import os, sys, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QC_SEC  = Path("${DATA_ROOT}/ST/results/step01_qc/section_h5ad")
COHORT  = Path("${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad")
OUT     = Path("${DATA_ROOT}/ST/results/step03_deconvolution/spatial_plots")
OUT.mkdir(parents=True, exist_ok=True)

PANELS = ["Malignant", "Macro_SPP1", "Neu_Inflammatory", "Neu_OSM_priming",
          "Neu_OSM_low", "Fibroblast", "Endothelial"]


def main():
    print(f"[load] {COHORT}")
    cohort = sc.read_h5ad(str(COHORT))
    abund = cohort.obsm["q05_cell_abundance"].copy()
    abund.index = cohort.obs_names
    print(f"   cohort: {cohort.shape}; q05 cols={abund.shape[1]}")

    samples = sorted(cohort.obs["sample"].unique().tolist())
    for s in samples:
        try:
            section_h5 = QC_SEC / f"{s}.h5ad"
            adata = sc.read_h5ad(str(section_h5))
            adata.var_names_make_unique()

            # cohort obs names are "<barcode>-<sample>"; map back to section barcodes
            mask = cohort.obs["sample"].values == s
            sub_idx = cohort.obs_names[mask]
            # strip the suffix "-<sample>"
            suffix = "-" + s
            section_bc = pd.Index([n[:-len(suffix)] if n.endswith(suffix) else n for n in sub_idx])

            sub = abund.loc[sub_idx].copy()
            sub.index = section_bc

            # Align to adata
            common = adata.obs_names.intersection(sub.index)
            print(f"[{s}] cohort spots={len(sub_idx)}, section spots={adata.n_obs}, common={len(common)}")
            adata = adata[common].copy()
            for col in PANELS:
                if col in sub.columns:
                    adata.obs[col] = sub.loc[adata.obs_names, col].values
                else:
                    adata.obs[col] = np.nan
                    print(f"   [warn] panel {col} missing in q05")

            # Plot 2x4 grid
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            axes = axes.flatten()
            for i, ct in enumerate(PANELS):
                ax = axes[i]
                try:
                    sc.pl.spatial(
                        adata, color=ct, library_id=s, ax=ax, show=False,
                        cmap="magma", size=1.4, frameon=False,
                        title=ct, colorbar_loc=None,
                    )
                except Exception as e:
                    ax.set_title(f"{ct} (err)")
                    ax.text(0.5, 0.5, f"{type(e).__name__}", ha="center", va="center")
            # last axis -> spot QC: total counts
            try:
                sc.pl.spatial(adata, color="total_counts", library_id=s, ax=axes[-1],
                              show=False, cmap="viridis", size=1.4, frameon=False,
                              title=f"total_counts (n_spots={adata.n_obs})", colorbar_loc=None)
            except Exception:
                axes[-1].axis("off")
            fig.suptitle(f"{s}  cell2location q05 abundance", fontsize=13)
            fig.tight_layout()
            out_png = OUT / f"{s}.png"
            fig.savefig(out_png, dpi=130, bbox_inches="tight")
            plt.close(fig)
            print(f"   saved {out_png}")
            del adata; gc.collect()
        except Exception as e:
            print(f"[ERROR] {s}: {type(e).__name__}: {e}\n{traceback.format_exc()}")

    print("[done] step03 plots")


if __name__ == "__main__":
    main()
