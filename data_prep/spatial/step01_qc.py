"""
Fig7 Step 1: Read E-MTAB-13530 Visium sections, run QC, filter LUAD tumor sections.

Inputs:
  ${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530/<sample>-filtered_feature_bc_matrix.h5
  ${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530/<sample>-spatial/{scalefactors,tissue_positions_list,tissue_hires_image,tissue_lowres_image}

Outputs:
  ${DATA_ROOT}/ST/results/step01_qc/section_h5ad/<sample>.h5ad   per-section QC'ed AnnData
  ${DATA_ROOT}/ST/results/step01_qc/luad_tumor_sections.h5ad     concatenated LUAD tumor cohort
  ${DATA_ROOT}/ST/results/step01_qc/qc_summary.csv               per-section pre/post QC stats
  ${DATA_ROOT}/ST/results/step01_qc/sample_metadata.csv          sdrf-derived per-section metadata

QC thresholds (typical Visium FF):
  - drop spots with total_counts < 500
  - drop spots with n_genes_by_counts < 200
  - drop spots with pct_counts_mt > 25
  - keep genes detected in >= 3 spots in the cohort
"""
from __future__ import annotations
import json, os, sys, gc
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from PIL import Image

ROOT = Path("${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530")
SDRF = ROOT / "E-MTAB-13530.sdrf.txt"
OUT  = Path("${DATA_ROOT}/ST/results/step01_qc")
OUT_H5AD = OUT / "section_h5ad"
OUT.mkdir(parents=True, exist_ok=True)
OUT_H5AD.mkdir(parents=True, exist_ok=True)

MIN_COUNTS = 500
MIN_GENES  = 200
MAX_MT_PCT = 25.0


def parse_metadata() -> pd.DataFrame:
    """Per-section metadata from sdrf, deduplicated to one row per Source Name."""
    df = pd.read_csv(SDRF, sep="\t", dtype=str)
    cols = {
        "Source Name": "section",
        "Characteristics[individual]": "patient",
        "Characteristics[disease]": "disease",
        "Characteristics[disease staging]": "stage",
        "Characteristics[sampling site]": "site",
        "Characteristics[sex]": "sex",
        "Characteristics[age]": "age",
    }
    sub = df[list(cols)].rename(columns=cols).drop_duplicates(subset=["section"])
    sub["is_tumor"]  = sub["site"].str.lower() == "tumor"
    sub["is_luad"]   = sub["disease"].str.lower() == "lung adenocarcinoma"
    sub["is_lusc"]   = sub["disease"].str.lower().str.contains("squamous", na=False)
    sub["is_normal"] = sub["disease"].str.lower() == "normal"
    sub["is_nsclc_unspec"] = sub["disease"].str.lower() == "non-small cell carcinoma"
    return sub.reset_index(drop=True)


def read_visium_section(sample: str) -> ad.AnnData:
    """Read flattened E-MTAB-13530 section into a Visium-compatible AnnData."""
    h5 = ROOT / f"{sample}-filtered_feature_bc_matrix.h5"
    sp_dir = ROOT / f"{sample}-spatial"
    adata = sc.read_10x_h5(str(h5))
    adata.var_names_make_unique()

    # Load scalefactors
    with open(sp_dir / "scalefactors_json.json") as f:
        sf = json.load(f)

    # Load tissue positions (no header in SpaceRanger 1.x output)
    tp = pd.read_csv(sp_dir / "tissue_positions_list.csv", header=None,
                     names=["barcode","in_tissue","array_row","array_col",
                            "pxl_row_in_fullres","pxl_col_in_fullres"])
    tp = tp.set_index("barcode")
    # restrict to spots actually in the matrix (already filtered_feature_bc_matrix is in-tissue)
    common = adata.obs_names.intersection(tp.index)
    adata = adata[common].copy()
    tp = tp.loc[adata.obs_names]
    adata.obs[["in_tissue","array_row","array_col"]] = tp[["in_tissue","array_row","array_col"]].values
    adata.obsm["spatial"] = tp[["pxl_col_in_fullres","pxl_row_in_fullres"]].to_numpy(dtype=float)

    # Load images
    hires_path = sp_dir / "tissue_hires_image.png"
    lowres_path = sp_dir / "tissue_lowres_image.png"
    images = {}
    if hires_path.exists():
        images["hires"] = np.asarray(Image.open(hires_path))
    if lowres_path.exists():
        images["lowres"] = np.asarray(Image.open(lowres_path))

    adata.uns["spatial"] = {
        sample: {
            "images": images,
            "scalefactors": sf,
            "metadata": {"chemistry_description": "Visium", "software_version": "spaceranger-1.1.0"},
        }
    }
    adata.obs["sample"] = sample
    return adata


def qc_one(adata: ad.AnnData, sample: str) -> tuple[ad.AnnData, dict]:
    pre_n_spots, pre_n_genes = adata.n_obs, adata.n_vars
    adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    keep = (
        (adata.obs["total_counts"]       >= MIN_COUNTS) &
        (adata.obs["n_genes_by_counts"]  >= MIN_GENES)  &
        (adata.obs["pct_counts_mt"]      <= MAX_MT_PCT)
    )
    adata = adata[keep].copy()
    summary = {
        "sample": sample,
        "pre_spots": pre_n_spots,
        "pre_genes": pre_n_genes,
        "post_spots": adata.n_obs,
        "post_genes": adata.n_vars,
        "median_counts": float(np.median(adata.obs["total_counts"])) if adata.n_obs else 0.0,
        "median_genes":  float(np.median(adata.obs["n_genes_by_counts"])) if adata.n_obs else 0.0,
        "median_mt_pct": float(np.median(adata.obs["pct_counts_mt"])) if adata.n_obs else 0.0,
    }
    return adata, summary


def main():
    meta = parse_metadata()
    meta.to_csv(OUT / "sample_metadata.csv", index=False)

    # Discover sections present on disk
    sections = sorted(p.name.replace("-filtered_feature_bc_matrix.h5","")
                      for p in ROOT.glob("*-filtered_feature_bc_matrix.h5"))
    meta = meta[meta["section"].isin(sections)].reset_index(drop=True)
    print(f"[info] {len(sections)} sections on disk; {len(meta)} matched in sdrf")

    # Process every section so we have full QC summary; only keep LUAD tumor in cohort
    summaries: list[dict] = []
    luad_tumor_h5ads: list[str] = []

    for _, row in meta.iterrows():
        sample = row["section"]
        try:
            print(f"[read] {sample}: {row['disease']} | {row['site']}")
            adata = read_visium_section(sample)
            adata, s = qc_one(adata, sample)
            s.update({k: row[k] for k in ["disease","site","patient","stage","sex","age"]})
            s["is_luad_tumor"] = bool(row["is_luad"] and row["is_tumor"])
            summaries.append(s)

            # Write every QCed section to its own h5ad
            out_h5 = OUT_H5AD / f"{sample}.h5ad"
            adata.write_h5ad(out_h5, compression="gzip")
            if s["is_luad_tumor"]:
                luad_tumor_h5ads.append(str(out_h5))
            del adata; gc.collect()
        except Exception as e:
            print(f"[error] {sample}: {type(e).__name__}: {e}", file=sys.stderr)
            summaries.append({"sample": sample, "error": f"{type(e).__name__}: {e}"})

    df_sum = pd.DataFrame(summaries)
    df_sum.to_csv(OUT / "qc_summary.csv", index=False)
    print(f"[info] qc_summary written: {OUT/'qc_summary.csv'}")

    # Build LUAD-tumor cohort (concatenate, batch_key='sample')
    if luad_tumor_h5ads:
        adatas = {Path(p).stem: sc.read_h5ad(p) for p in luad_tumor_h5ads}
        # Drop heavy uns['spatial'] images for the cohort file (per-sample h5ads keep them)
        for k, a in adatas.items():
            if "spatial" in a.uns:
                a.uns["spatial"] = {k: {"scalefactors": a.uns["spatial"][k]["scalefactors"]}}
        cohort = ad.concat(adatas, label="sample", join="outer", index_unique="-", merge="unique")
        cohort.obs_names_make_unique()
        cohort_path = OUT / "luad_tumor_sections.h5ad"
        cohort.write_h5ad(cohort_path, compression="gzip")
        print(f"[info] cohort written: {cohort_path}  ({cohort.n_obs} spots × {cohort.n_vars} genes)")

    # Print headline report
    df_lt = df_sum[df_sum.get("is_luad_tumor", False) == True]
    print("\n=== LUAD tumor sections (post-QC) ===")
    if len(df_lt):
        print(df_lt[["sample","patient","post_spots","post_genes","median_counts","median_genes","median_mt_pct"]].to_string(index=False))
        print(f"\ntotal LUAD tumor sections: {len(df_lt)}")
        print(f"total LUAD tumor spots:    {int(df_lt['post_spots'].sum())}")
    else:
        print("no LUAD tumor sections found!")


if __name__ == "__main__":
    main()
