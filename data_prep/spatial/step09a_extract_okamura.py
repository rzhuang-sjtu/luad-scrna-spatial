"""
Step 9a: extract LUAD No.1-5 FF Visium tarballs, then QC + build joint AnnData.

Inputs:
  ${DATA_ROOT}/ST/Okamura 2024/Visium_FF_LUAD_No_{1,2,3,4,5}.tar.gz
Outputs:
  ${DATA_ROOT}/ST/results/step09_okamura_validation/raw/LUAD_No_<N>/   extracted SpaceRanger outs
  ${DATA_ROOT}/ST/results/step09_okamura_validation/section_h5ad/<sample>.h5ad
  ${DATA_ROOT}/ST/results/step09_okamura_validation/cohort.h5ad         joint cohort (ready for c2l)
  ${DATA_ROOT}/ST/results/step09_okamura_validation/qc_summary.csv
  ${DATA_ROOT}/ST/results/step09_okamura_validation/memory_estimate.txt

QC: same thresholds as Step 1 (min_counts>=500, min_genes>=200, mt%<=25).
"""
from __future__ import annotations
import os, sys, time, gc, json, tarfile
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from PIL import Image

OKAMURA = Path("${DATA_ROOT}/ST/Okamura 2024")
# 16 invasive LUAD Visium sections from 8 patients (Takano 2024).
# 8 fresh-frozen (FF) + 8 formalin-fixed paraffin-embedded (FFPE) multi-section
# from patients No.2, No.3, No.4. Sample IDs are unique across FF/FFPE.
#
# All 16 are unpacked here, but only the 8 FF sections enter the analysis: the
# validation cohort in the paper is FF-only, and the downstream steps filter on
# that list (see TAKANO_FF in analysis/supplementary_tables/). The FFPE sections
# were unpacked while deciding whether to include them and are kept so the
# decision is visible rather than silent.
SAMPLES = [
    # FF
    ("LUAD_No_1",  "Visium_FF_LUAD_No_1.tar.gz",   "FF"),
    ("LUAD_No_2",  "Visium_FF_LUAD_No_2.tar.gz",   "FF"),
    ("LUAD_No_3",  "Visium_FF_LUAD_No_3.tar.gz",   "FF"),
    ("LUAD_No_4",  "Visium_FF_LUAD_No_4.tar.gz",   "FF"),
    ("LUAD_No_5",  "Visium_FF_LUAD_No_5.tar.gz",   "FF"),
    ("LUAD_No_14", "Visium_FF_LUAD_No_14.tar.gz",  "FF"),
    ("LUAD_No_16", "Visium_FF_LUAD_No_16.tar.gz",  "FF"),
    ("LUAD_No_17", "Visium_FF_LUAD_No_17.tar.gz",  "FF"),
    # FFPE — multi-section from same patient
    ("LUAD_No_2A", "Visium_FFPE_LUAD_No_2A.tar.gz", "FFPE"),
    ("LUAD_No_2B", "Visium_FFPE_LUAD_No_2B.tar.gz", "FFPE"),
    ("LUAD_No_2C", "Visium_FFPE_LUAD_No_2C.tar.gz", "FFPE"),
    ("LUAD_No_2D", "Visium_FFPE_LUAD_No_2D.tar.gz", "FFPE"),
    ("LUAD_No_3A", "Visium_FFPE_LUAD_No_3A.tar.gz", "FFPE"),
    ("LUAD_No_3B", "Visium_FFPE_LUAD_No_3B.tar.gz", "FFPE"),
    ("LUAD_No_4C", "Visium_FFPE_LUAD_No_4C.tar.gz", "FFPE"),
    ("LUAD_No_4D", "Visium_FFPE_LUAD_No_4D.tar.gz", "FFPE"),
]
ROOT    = Path("${DATA_ROOT}/ST/results/step09_okamura_validation")
RAW     = ROOT / "raw"
SEC_OUT = ROOT / "section_h5ad"
for d in (ROOT, RAW, SEC_OUT):
    d.mkdir(parents=True, exist_ok=True)
LOG = ROOT / "run.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

MIN_COUNTS, MIN_GENES, MAX_MT = 500, 200, 25.0


def extract_one(sample: str, tar_name: str):
    tar = OKAMURA / tar_name
    sd = RAW / sample
    h5 = sd / "filtered_feature_bc_matrix.h5"
    sp = sd / "spatial"
    if h5.exists() and sp.exists():
        log(f"  {sample}: already extracted")
        return sd
    if not tar.exists():
        raise FileNotFoundError(f"{tar} not found")
    log(f"  {sample}: extracting {tar.name} -> {sd.parent}")
    sd.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar, "r:gz") as t:
        # The tarballs use sample-specific top-level dirs (e.g.
        # "Visium_FFPE_LUAD_No_2A/outs/..."). Strip the leading dir component
        # so we always end up with sd/{filtered_feature_bc_matrix.h5,spatial/}
        wanted = []
        for m in t.getmembers():
            n = m.name
            # keep raw matrix + spatial dir + metrics summary
            keep = (n.endswith("filtered_feature_bc_matrix.h5") or
                    "/spatial/" in n or n.endswith("/spatial") or
                    n.endswith("metrics_summary.csv"))
            if not keep:
                continue
            # Drop top-level dir, also drop "outs/" if present
            parts = n.split("/")
            if parts and parts[0] != sample:
                # strip arbitrary top-level dir
                parts = parts[1:]
            if parts and parts[0] == "outs":
                parts = parts[1:]
            if not parts:
                continue
            m.name = "/".join([sample] + parts)
            wanted.append(m)
        for m in wanted:
            try:
                t.extract(m, path=sd.parent, filter="data")
            except TypeError:
                t.extract(m, path=sd.parent)
    log(f"  {sample}: extracted {len(wanted)} entries")
    return sd


def read_visium_section(sample: str, sample_dir: Path) -> ad.AnnData:
    h5 = sample_dir / "filtered_feature_bc_matrix.h5"
    sp_dir = sample_dir / "spatial"
    adata = sc.read_10x_h5(str(h5))
    adata.var_names_make_unique()
    with open(sp_dir / "scalefactors_json.json") as f:
        sf = json.load(f)
    tp_path = sp_dir / "tissue_positions_list.csv"
    if not tp_path.exists():
        tp_path = sp_dir / "tissue_positions.csv"
    # No header in SpaceRanger 1.x (positions_list.csv); has header in 2.x (positions.csv)
    if "tissue_positions_list" in tp_path.name:
        tp = pd.read_csv(tp_path, header=None,
                         names=["barcode","in_tissue","array_row","array_col",
                                "pxl_row_in_fullres","pxl_col_in_fullres"])
    else:
        tp = pd.read_csv(tp_path)
        tp.columns = ["barcode","in_tissue","array_row","array_col",
                      "pxl_row_in_fullres","pxl_col_in_fullres"][: len(tp.columns)]
    tp = tp.set_index("barcode")
    common = adata.obs_names.intersection(tp.index)
    adata = adata[common].copy()
    tp = tp.loc[adata.obs_names]
    adata.obs[["in_tissue","array_row","array_col"]] = tp[["in_tissue","array_row","array_col"]].values
    adata.obsm["spatial"] = tp[["pxl_col_in_fullres","pxl_row_in_fullres"]].to_numpy(dtype=float)
    images = {}
    for k, fn in (("hires","tissue_hires_image.png"), ("lowres","tissue_lowres_image.png")):
        p = sp_dir / fn
        if p.exists(): images[k] = np.asarray(Image.open(p))
    adata.uns["spatial"] = {sample: {"images": images, "scalefactors": sf,
                                     "metadata": {"chemistry_description": "Visium",
                                                  "software_version": "spaceranger"}}}
    adata.obs["sample"] = sample
    return adata


def qc_one(adata: ad.AnnData, sample: str) -> tuple[ad.AnnData, dict]:
    pre_n, pre_g = adata.n_obs, adata.n_vars
    adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    keep = ((adata.obs["total_counts"] >= MIN_COUNTS) &
            (adata.obs["n_genes_by_counts"] >= MIN_GENES) &
            (adata.obs["pct_counts_mt"] <= MAX_MT))
    adata = adata[keep].copy()
    return adata, {
        "sample": sample, "pre_spots": pre_n, "pre_genes": pre_g,
        "post_spots": adata.n_obs, "post_genes": adata.n_vars,
        "median_counts": float(np.median(adata.obs["total_counts"])) if adata.n_obs else 0.0,
        "median_genes":  float(np.median(adata.obs["n_genes_by_counts"])) if adata.n_obs else 0.0,
        "median_mt_pct": float(np.median(adata.obs["pct_counts_mt"])) if adata.n_obs else 0.0,
    }


def main():
    log(f"samples: {[s[0] for s in SAMPLES]}")
    summaries = []
    adatas = {}
    for s, tar_name, prep in SAMPLES:
        try:
            sd = extract_one(s, tar_name)
        except FileNotFoundError as e:
            log(f"  [SKIP] {s}: {e}")
            continue
        adata = read_visium_section(s, sd)
        adata.obs["preparation"] = prep   # FF / FFPE
        adata, sm = qc_one(adata, s)
        sm["preparation"] = prep
        log(f"  {s} ({prep}): post-QC {adata.shape}; median_counts={sm['median_counts']:.0f}; median_genes={sm['median_genes']:.0f}")
        adata.write_h5ad(str(SEC_OUT / f"{s}.h5ad"), compression="gzip")
        # drop heavy uns for cohort
        if "spatial" in adata.uns and isinstance(adata.uns["spatial"], dict):
            for k, v in list(adata.uns["spatial"].items()):
                adata.uns["spatial"][k] = {kk: v[kk] for kk in ("scalefactors",) if kk in v}
        adatas[s] = adata
        summaries.append(sm)

    df = pd.DataFrame(summaries)
    df.to_csv(ROOT / "qc_summary.csv", index=False)
    log(f"\n=== qc summary ===\n{df.to_string(index=False)}")

    cohort = ad.concat(adatas, label="sample", join="outer", index_unique="-",
                       merge="unique", uns_merge="unique")
    cohort.obs_names_make_unique()
    cohort.write_h5ad(str(ROOT / "cohort.h5ad"), compression="gzip")
    log(f"cohort: {cohort.n_obs} spots × {cohort.n_vars} genes  -> cohort.h5ad")

    # Memory estimate for cell2location
    n_spots = cohort.n_obs
    inf_aver = pd.read_csv("${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv", index_col=0)
    n_genes_overlap = len(set(cohort.var_names) & set(inf_aver.index))
    n_cell_types = inf_aver.shape[1]
    # Rough estimate: c2l full-batch holds ~3 tensors of shape (n_spots, n_cell_types) +
    # (n_spots, n_genes) + autograd doubles → ~4× peak.
    #   spot×ct ≈ n_spots * n_cell_types * 4  (float32)
    #   spot×gene ≈ n_spots * n_genes * 4
    #   ct×gene  ≈ n_cell_types * n_genes * 4
    # Total ≈ 4 * (spot×ct + spot×gene + ct×gene) ; multiplied by ~6 for autograd/posterior.
    spot_ct  = n_spots * n_cell_types * 4
    spot_g   = n_spots * n_genes_overlap * 4
    ct_g     = n_cell_types * n_genes_overlap * 4
    rough_gb = (spot_ct + spot_g + ct_g) * 6 / 1e9
    msg = (
        f"n_spots={n_spots}\n"
        f"n_genes_overlap_with_ref={n_genes_overlap}\n"
        f"n_cell_types={n_cell_types}\n"
        f"rough_peak_GPU_GB_estimate={rough_gb:.2f}\n"
        f"recommend_full_batch={'yes' if rough_gb < 14 else 'no, fall back batch_size=2500'}\n"
    )
    open(ROOT / "memory_estimate.txt", "w").write(msg)
    log("\n" + msg)


if __name__ == "__main__":
    main()
