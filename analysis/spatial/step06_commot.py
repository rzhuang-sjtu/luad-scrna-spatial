"""
Step 6: COMMOT spatial communication on OSM and IL1 pathways, per section.

For each LUAD tumor section:
  - load QC'ed Visium AnnData (per-section)
  - normalize log1p CPM
  - run commot.tl.spatial_communication with CellChat human DB filtered to OSM + IL1 pathways
  - run commot.tl.communication_direction to get vector field
  - save section h5ad with communication results in obsm / obsp
  - plot 2x3 panel: OSM/IL1 sender + receiver + total scores

Outputs:
  ${DATA_ROOT}/ST/results/step06_commot/section_h5ad/<sample>.h5ad
  ${DATA_ROOT}/ST/results/step06_commot/spatial_plots/<sample>.png
  ${DATA_ROOT}/ST/results/step06_commot/per_sample_pathway_summary.csv
"""
from __future__ import annotations
import os, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QC_DIR = Path("${DATA_ROOT}/ST/results/step01_qc")
QC_SEC = QC_DIR / "section_h5ad"
OUT    = Path("${DATA_ROOT}/ST/results/step06_commot")
SEC_OUT = OUT / "section_h5ad"
PLOTS  = OUT / "spatial_plots"
for d in (OUT, SEC_OUT, PLOTS):
    d.mkdir(parents=True, exist_ok=True)
LOG = OUT / "run.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

PATHWAYS = ["OSM", "IL1"]   # target pathways for Fig 7 E-G
DIS_THR = 500.0            # pixel-space neighborhood radius for Visium hi-res coords


def main():
    qc = pd.read_csv(QC_DIR / "qc_summary.csv")
    samples = sorted(qc[qc["is_luad_tumor"] == True]["sample"].tolist())
    log(f"sections: {len(samples)} -> {samples}")

    import commot as ct
    log(f"commot version: {ct.__version__ if hasattr(ct,'__version__') else 'unknown'}")

    # Load CellChat human LR database
    df_lr_full = ct.pp.ligand_receptor_database(species="human", signaling_type=None, database="CellChat")
    log(f"cellchat human LR db: {df_lr_full.shape}; example rows:\n{df_lr_full.head(3).to_string()}")
    # CellChat DB columns: [0]=ligand, [1]=receptor, [2]=pathway, [3]=signaling_type
    df_lr = df_lr_full[df_lr_full.iloc[:, 2].isin(PATHWAYS)].copy()
    log(f"filtered to {PATHWAYS}: {df_lr.shape}\n{df_lr.to_string()}")
    # Save for record
    df_lr.to_csv(OUT / "lr_pairs_used.csv", index=False)

    summary_rows = []
    for i, s in enumerate(samples, 1):
        try:
            t0 = time.time()
            log(f"[{i}/{len(samples)}] {s}: loading section ...")
            sec = sc.read_h5ad(str(QC_SEC / f"{s}.h5ad"))
            sec.var_names_make_unique()
            sec.layers["counts"] = sec.X.copy()
            sc.pp.normalize_total(sec, target_sum=1e4)
            sc.pp.log1p(sec)

            # Run spatial communication
            log(f"    spatial_communication: dis_thr={DIS_THR}, heteromeric=True ...")
            ct.tl.spatial_communication(
                sec, database_name="cellchat", df_ligrec=df_lr,
                dis_thr=DIS_THR, heteromeric=True, pathway_sum=True,
            )

            # Direction vector fields
            for p in PATHWAYS:
                try:
                    ct.tl.communication_direction(sec, database_name="cellchat",
                                                   pathway_name=p, k=5)
                except Exception as e:
                    log(f"    [warn] communication_direction({p}): {type(e).__name__}: {e}")

            # Sanity: capture available obsm / obs keys
            commot_obs   = [c for c in sec.obs.columns if c.startswith("s-cellchat") or c.startswith("r-cellchat")]
            commot_obsm  = [k for k in sec.obsm.keys()  if "cellchat" in k.lower() or "commot" in k.lower()]
            commot_obsp  = [k for k in sec.obsp.keys()  if "cellchat" in k.lower() or "commot" in k.lower()]
            log(f"    commot obs cols: {len(commot_obs)}; obsm keys: {len(commot_obsm)}; obsp keys: {len(commot_obsp)}")

            # Save
            sec.write_h5ad(str(SEC_OUT / f"{s}.h5ad"), compression="gzip")

            # Plot 2x3: per pathway sender/receiver + total
            fig, axes = plt.subplots(2, 3, figsize=(15, 9))
            for r, p in enumerate(PATHWAYS):
                # try potential commot obs keys
                candidates = {
                    f"s-cellchat-{p}": axes[r, 0],   # sender (outgoing) — common naming
                    f"r-cellchat-{p}": axes[r, 1],
                    f"s-cellchat-pathway_{p}": axes[r, 0],
                    f"r-cellchat-pathway_{p}": axes[r, 1],
                }
                # more flexible: detect actual columns
                send_col = next((c for c in sec.obs.columns if c.startswith("s-cellchat") and p in c), None)
                recv_col = next((c for c in sec.obs.columns if c.startswith("r-cellchat") and p in c), None)
                tot_col  = None
                if send_col and recv_col:
                    sec.obs[f"total-{p}"] = sec.obs[send_col].fillna(0) + sec.obs[recv_col].fillna(0)
                    tot_col = f"total-{p}"
                ax = axes[r, 0]
                if send_col:
                    sc.pl.spatial(sec, color=send_col, library_id=s, ax=ax, show=False,
                                  cmap="magma", size=1.4, frameon=False, title=f"{p} sender", colorbar_loc="right")
                else:
                    ax.set_title(f"{p} sender (NA)"); ax.axis("off")
                ax = axes[r, 1]
                if recv_col:
                    sc.pl.spatial(sec, color=recv_col, library_id=s, ax=ax, show=False,
                                  cmap="magma", size=1.4, frameon=False, title=f"{p} receiver", colorbar_loc="right")
                else:
                    ax.set_title(f"{p} receiver (NA)"); ax.axis("off")
                ax = axes[r, 2]
                if tot_col:
                    sc.pl.spatial(sec, color=tot_col, library_id=s, ax=ax, show=False,
                                  cmap="magma", size=1.4, frameon=False, title=f"{p} total", colorbar_loc="right")
                else:
                    ax.set_title(f"{p} total (NA)"); ax.axis("off")
            fig.suptitle(f"{s}  COMMOT (CellChat) — {', '.join(PATHWAYS)}", fontsize=12)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}.png", dpi=130, bbox_inches="tight")
            plt.close(fig)

            # Summary
            row = {"sample": s, "n_spots": sec.n_obs, "elapsed_min": (time.time()-t0)/60}
            for p in PATHWAYS:
                send_col = next((c for c in sec.obs.columns if c.startswith("s-cellchat") and p in c), None)
                recv_col = next((c for c in sec.obs.columns if c.startswith("r-cellchat") and p in c), None)
                row[f"{p}_send_mean"] = float(sec.obs[send_col].mean()) if send_col else np.nan
                row[f"{p}_recv_mean"] = float(sec.obs[recv_col].mean()) if recv_col else np.nan
            summary_rows.append(row)
            log(f"    {s} done in {row['elapsed_min']:.1f} min")
            del sec; gc.collect()
        except Exception as e:
            log(f"[ERROR] {s}: {type(e).__name__}: {e}\n{traceback.format_exc()}")
            summary_rows.append({"sample": s, "error": f"{type(e).__name__}: {e}"})

    pd.DataFrame(summary_rows).to_csv(OUT / "per_sample_pathway_summary.csv", index=False)
    log("[done]")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
