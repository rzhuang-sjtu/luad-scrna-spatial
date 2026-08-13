"""
Step A: export per-section R-friendly CSVs for both cohorts.

For each section in E-MTAB-13530 (12) and Okamura (8):
  per_section/<sample>.csv
    spot_id, spatial1, spatial2, in_tissue, sample, cohort,
    [28 cell type abundances]   prefixed: ct_<celltype>
    [MP1-5 score]                  MP1_score..MP5_score
    dominant_MP_4
    [14 PROGENy pathway scores]    progeny_<pathway>
    neu_total
    z_nfkb, z_neu, roi
    [gene expression for FIG7 markers, log-normalized]
        gex_OSM, gex_IL1B, gex_ATF3, gex_FOSB, gex_JUN, gex_JUNB, gex_NFKBIA,
        gex_KRT7, gex_LAMC2, gex_TGFB1, gex_FOS, gex_HIF1A, gex_IL1A, gex_IL1R1
    [COMMOT pathway-level scores]  commot_s_OSM, commot_r_OSM, commot_total_OSM,
                                   commot_s_IL1, commot_r_IL1, commot_total_IL1

Plus shared tables:
  commot_per_sample_summary.csv  (cohort, sample, OSM/IL1 send/recv mean)
  misty_aggregated_importance.csv (cohort, view, target, predictor, mean_importance)
  roi_vs_nonroi_aggregate.csv     (cohort, metric, type, mean_roi, mean_nonroi, delta)
  qc_section_metadata.csv         (cohort, sample, n_spots, ...)
"""
from __future__ import annotations
import os, time, json, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp

OUT = Path("${DATA_ROOT}/ST/results/r_data")
OUT_PERSECTION = OUT / "per_section"
OUT.mkdir(parents=True, exist_ok=True)
OUT_PERSECTION.mkdir(parents=True, exist_ok=True)
LOG = OUT / "export.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

GENES = ["OSM","IL1B","IL1A","IL1R1","OSMR","LIFR","ATF3","FOSB","JUN","JUNB","FOS",
         "NFKBIA","KRT7","LAMC2","TGFB1","HIF1A","SPP1","CD44","CXCL8","CXCL1"]

EMTAB_COHORT = Path("${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad")
EMTAB_COMMOT_SEC = Path("${DATA_ROOT}/ST/results/step06_commot/section_h5ad")

OKA_COHORT = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/cohort_with_roi.h5ad")
OKA_COMMOT_SEC = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/commot")

EMTAB_MISTY = Path("${DATA_ROOT}/ST/results/step07_misty/aggregated_importance.csv")
OKA_MISTY   = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/misty/aggregated_importance.csv")
EMTAB_COMMOT_SUM = Path("${DATA_ROOT}/ST/results/step06_commot/per_sample_pathway_summary.csv")
OKA_COMMOT_SUM   = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/per_sample_pathway_summary.csv")
EMTAB_ROI_AGG = Path("${DATA_ROOT}/ST/results/step08_roi/roi_vs_nonroi_aggregate.csv")
OKA_ROI_AGG   = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/roi_vs_nonroi_aggregate.csv")
EMTAB_QC = Path("${DATA_ROOT}/ST/results/step01_qc/qc_summary.csv")
OKA_QC   = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/qc_summary.csv")


def export_cohort(cohort_path: Path, commot_dir: Path, cohort_label: str):
    log(f"loading {cohort_path}")
    a = sc.read_h5ad(str(cohort_path))
    log(f"   shape={a.shape}; cohort={cohort_label}")
    abund = a.obsm["q05_cell_abundance"].copy()
    abund.index = a.obs_names
    samples = sorted(a.obs["sample"].unique().tolist())

    # gene expression columns (log-normalized in cohort.X already)
    genes_present = [g for g in GENES if g in a.var_names]
    log(f"   genes present: {genes_present}")
    X = a[:, genes_present].X
    if sp.issparse(X): X = X.toarray()
    gex = pd.DataFrame(X, index=a.obs_names, columns=[f"gex_{g}" for g in genes_present])

    # PROGENy obs cols
    progeny_cols = [c for c in a.obs.columns if c.startswith("progeny_")]
    mp_cols = [c for c in a.obs.columns if c.startswith("MP") and c.endswith("_score")]

    for s in samples:
        mask = a.obs["sample"].values == s
        sub_idx = a.obs_names[mask]

        df = pd.DataFrame(index=sub_idx)
        df.index.name = "spot_id"
        # spatial coords
        coords = a.obsm["spatial"][mask]
        df["spatial1"] = coords[:, 0]
        df["spatial2"] = coords[:, 1]
        df["sample"]   = s
        df["cohort"]   = cohort_label
        df["in_tissue"] = a.obs.loc[sub_idx, "in_tissue"].values if "in_tissue" in a.obs.columns else 1
        # cell type abundances (28 cols)
        ab_sub = abund.loc[sub_idx].copy()
        ab_sub.columns = [f"ct_{c}" for c in ab_sub.columns]
        df = df.join(ab_sub)
        # MP scores
        for col in mp_cols:
            df[col] = a.obs.loc[sub_idx, col].values
        if "dominant_MP_4" in a.obs.columns:
            df["dominant_MP_4"] = a.obs.loc[sub_idx, "dominant_MP_4"].astype(str).values
        # PROGENy
        for col in progeny_cols:
            df[col] = a.obs.loc[sub_idx, col].values
        # neu_total / z / roi
        for col in ["neu_total","z_nfkb","z_neu","roi"]:
            if col in a.obs.columns:
                df[col] = a.obs.loc[sub_idx, col].values
        if "roi" in df.columns:
            df["roi"] = df["roi"].astype(int)
        # gene expression
        df = df.join(gex.loc[sub_idx])

        # COMMOT obs columns from per-section h5ad
        commot_h5 = commot_dir / f"{s}.h5ad"
        if commot_h5.exists():
            csec = sc.read_h5ad(str(commot_h5), backed="r")
            commot_cols = ["s_OSM","r_OSM","total_OSM","s_IL1","r_IL1","total_IL1"]
            obs_cols_present = [c for c in commot_cols if c in csec.obs.columns]
            if obs_cols_present:
                # commot uses raw barcodes (no -<sample> suffix)
                comm_obs = csec.obs[obs_cols_present].copy()
                # cohort obs_names look like "<barcode>-<sample>"
                strip_suffix = "-" + s
                section_bc = [n[:-len(strip_suffix)] if n.endswith(strip_suffix) else n for n in df.index]
                # remap
                df_remap = pd.DataFrame(index=section_bc)
                df_remap.index.name = "section_barcode"
                for c in obs_cols_present:
                    df_remap[f"commot_{c}"] = comm_obs.reindex(df_remap.index)[c].values
                df_remap.index = df.index
                df = df.join(df_remap)
            csec.file.close()

        out_csv = OUT_PERSECTION / f"{cohort_label}__{s}.csv"
        df.to_csv(out_csv)
        log(f"  {cohort_label}/{s}: {df.shape} -> {out_csv.name}")
    a.file.close() if hasattr(a, "file") else None


def main():
    log("=== Step A: export r_data ===")

    # E-MTAB
    export_cohort(EMTAB_COHORT, EMTAB_COMMOT_SEC, "EMTAB13530")
    # Okamura
    export_cohort(OKA_COHORT,  OKA_COMMOT_SEC,  "Okamura")

    # Combined helper tables
    log("combining MISTy aggregated importance ...")
    misty_emtab = pd.read_csv(EMTAB_MISTY); misty_emtab["cohort"] = "EMTAB13530"
    misty_oka   = pd.read_csv(OKA_MISTY);   misty_oka["cohort"]   = "Okamura"
    misty = pd.concat([misty_emtab, misty_oka], ignore_index=True)
    misty.to_csv(OUT / "misty_aggregated_importance.csv", index=False)
    log(f"   {misty.shape}")

    log("combining COMMOT per-sample summary ...")
    co_emtab = pd.read_csv(EMTAB_COMMOT_SUM); co_emtab["cohort"] = "EMTAB13530"
    co_oka   = pd.read_csv(OKA_COMMOT_SUM);   co_oka["cohort"]   = "Okamura"
    co = pd.concat([co_emtab, co_oka], ignore_index=True)
    co.to_csv(OUT / "commot_per_sample_summary.csv", index=False)
    log(f"   {co.shape}")

    log("combining ROI aggregate ...")
    r_emtab = pd.read_csv(EMTAB_ROI_AGG); r_emtab["cohort"] = "EMTAB13530"
    r_oka   = pd.read_csv(OKA_ROI_AGG);   r_oka["cohort"]   = "Okamura"
    r = pd.concat([r_emtab, r_oka], ignore_index=True)
    r.to_csv(OUT / "roi_vs_nonroi_aggregate.csv", index=False)
    log(f"   {r.shape}")

    log("combining QC metadata ...")
    qc_emtab = pd.read_csv(EMTAB_QC); qc_emtab["cohort"] = "EMTAB13530"
    qc_oka   = pd.read_csv(OKA_QC);   qc_oka["cohort"]   = "Okamura"
    qc = pd.concat([qc_emtab, qc_oka], ignore_index=True, sort=False)
    qc.to_csv(OUT / "qc_section_metadata.csv", index=False)
    log(f"   {qc.shape}")

    # Section list with key recommendation flags
    log("computing recommended section per panel ...")
    # for each cohort, find best Macro_SPP1 + best COMMOT + most ROI
    rec = []
    for cohort_label, qc_csv, commot_sum, roi_path in [
        ("EMTAB13530", EMTAB_QC, EMTAB_COMMOT_SUM, "${DATA_ROOT}/ST/results/step08_roi/roi_summary.csv"),
        ("Okamura",   OKA_QC,    OKA_COMMOT_SUM,  "${DATA_ROOT}/ST/results/step09_okamura_validation/roi_summary.csv"),
    ]:
        co = pd.read_csv(commot_sum)
        roi = pd.read_csv(roi_path)
        merged = co.merge(roi, on="sample", how="outer")
        merged["cohort"] = cohort_label
        rec.append(merged)
    rec_df = pd.concat(rec, ignore_index=True)
    rec_df.to_csv(OUT / "panel_selection_helper.csv", index=False)
    log(f"   panel_selection_helper.csv -> {rec_df.shape}")
    print(rec_df.to_string())

    log("[done] r_data export")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
