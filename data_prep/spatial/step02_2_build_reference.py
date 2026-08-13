"""
Step 2.2: Build a unified scRNA-seq reference for cell2location (Fig 7).

Strategy:
  1. Load luad_merged_raw.h5ad (raw counts, 853469 × 9881)
  2. Pull obs labels (celltype_coarse, tissue_type) from luad_merged_annotated.h5ad
  3. Filter to lung-relevant tissues (Primary_Tumor + Adjacent_Normal + Normal_Lung)
  4. Initialize cell_type_fine from celltype_coarse (8 broad classes)
  5. Override Myeloid cells with myeloid_subtype from luad_myeloid.h5ad
  6. Within Macro_general, derive Macro_SPP1 (sc_Macro_SPP1 ≥ Q3) and Macro_FOLR2 (sc_Macro_FOLR2 ≥ Q3);
     remaining Macro_general cells stay labeled Macro_general
  7. Override neutrophils with neu_subtype from luad_neutrophil_own_annotated.h5ad
     (drop Neu_unclassified)
  8. Override malignant Epithelial cells with "Malignant" from luad_malignant_scored.h5ad (where malignant=True)
  9. Drop cell types with <50 cells; cap each cell type at MAX_PER_TYPE for tractable c2l training
  10. Save: unified_reference.h5ad (raw counts, .obs.cell_type_fine), cell_counts_report.csv

Caveats:
  - merged_raw .X is csr_matrix float32 with integer values (verified in step02_0)
  - obs_names should be consistent across all h5ads (they all derive from the same merged source)
"""
from __future__ import annotations
import os, gc
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import anndata as ad

PROC = Path(os.path.expanduser("~/luad/data/processed"))
OUT  = Path("${DATA_ROOT}/ST/results/step02_reference")
OUT.mkdir(parents=True, exist_ok=True)

MAX_PER_TYPE = 15000           # cap per cell type to keep c2l NB regression tractable
MIN_PER_TYPE = 50              # drop tiny categories
KEEP_TISSUES = {"Primary_Tumor", "Adjacent_Normal", "Normal_Lung"}

RNG = np.random.default_rng(0)


def load_obs_only(path: Path) -> pd.DataFrame:
    a = sc.read_h5ad(str(path), backed="r")
    obs = a.obs.copy()
    obs.index = obs.index.astype(str)
    a.file.close()
    return obs


def main():
    # 1. Load raw counts (full)
    print("[1] loading luad_merged_raw.h5ad ...")
    raw = sc.read_h5ad(str(PROC / "luad_merged_raw.h5ad"))
    print(f"   raw: {raw.shape}, X={type(raw.X).__name__} dtype={raw.X.dtype}")
    raw.var_names_make_unique()

    # 2. Pull broad annotation from annotated h5ad (same shape, same obs order)
    print("[2] loading celltype_coarse from luad_merged_annotated.h5ad ...")
    anno = sc.read_h5ad(str(PROC / "luad_merged_annotated.h5ad"), backed="r")
    assert anno.n_obs == raw.n_obs, "annotated/raw obs mismatch"
    anno_obs = anno.obs[["celltype_coarse", "celltype_marker", "tissue_type",
                         "patient_id", "sample_id", "dataset", "tissue_type_original"]].copy()
    anno.file.close()
    # join by index — confirm alignment first
    same = (raw.obs_names == anno_obs.index).mean()
    print(f"   raw vs annotated obs_names alignment: {same*100:.2f}%")
    if same < 1.0:
        # reorder
        anno_obs = anno_obs.reindex(raw.obs_names)
    for c in anno_obs.columns:
        raw.obs[c] = anno_obs[c].values

    # 3. Filter to lung-relevant tissues
    print(f"[3] filtering by tissue_type ∈ {KEEP_TISSUES} ...")
    keep_mask = raw.obs["tissue_type"].isin(KEEP_TISSUES).values
    print(f"   keep {keep_mask.sum()}/{raw.n_obs} cells")
    raw = raw[keep_mask].copy()
    gc.collect()

    # 4. Initialize cell_type_fine from celltype_coarse
    raw.obs["cell_type_fine"] = raw.obs["celltype_coarse"].astype(str)

    # 5. Pull myeloid_subtype from luad_myeloid.h5ad
    print("[5] merging myeloid_subtype ...")
    myo = load_obs_only(PROC / "luad_myeloid.h5ad")
    cols_needed = ["myeloid_subtype", "sc_Macro_SPP1", "sc_Macro_FOLR2"]
    myo = myo[cols_needed]
    common = raw.obs_names.intersection(myo.index)
    print(f"   myeloid cells overlap: {len(common)}")
    raw.obs.loc[common, "cell_type_fine"] = myo.loc[common, "myeloid_subtype"].astype(str).values
    # add scores aligned to raw
    raw.obs["sc_Macro_SPP1"]  = pd.Series(np.nan, index=raw.obs_names, dtype=float)
    raw.obs["sc_Macro_FOLR2"] = pd.Series(np.nan, index=raw.obs_names, dtype=float)
    raw.obs.loc[common, "sc_Macro_SPP1"]  = myo.loc[common, "sc_Macro_SPP1"].values
    raw.obs.loc[common, "sc_Macro_FOLR2"] = myo.loc[common, "sc_Macro_FOLR2"].values

    # 6. Derive Macro_SPP1 and Macro_FOLR2 inside Macro_general
    print("[6] deriving Macro_SPP1 and Macro_FOLR2 from Macro_general by score Q3 ...")
    is_general = (raw.obs["cell_type_fine"].values == "Macro_general")
    if is_general.sum():
        spp1_q3  = np.nanquantile(raw.obs.loc[is_general, "sc_Macro_SPP1"], 0.75)
        folr2_q3 = np.nanquantile(raw.obs.loc[is_general, "sc_Macro_FOLR2"], 0.75)
        print(f"   Macro_general n={int(is_general.sum())}, SPP1 Q3={spp1_q3:.3f}, FOLR2 Q3={folr2_q3:.3f}")
        # Choose the higher-scoring one if both pass threshold (avoid double-assign)
        sub = raw.obs.loc[is_general, ["sc_Macro_SPP1", "sc_Macro_FOLR2"]].copy()
        spp1_pass  = sub["sc_Macro_SPP1"]  >= spp1_q3
        folr2_pass = sub["sc_Macro_FOLR2"] >= folr2_q3
        # priority: whichever score is larger (z-normalized for fairness)
        spp1_z  = (sub["sc_Macro_SPP1"]  - sub["sc_Macro_SPP1"].mean())  / sub["sc_Macro_SPP1"].std()
        folr2_z = (sub["sc_Macro_FOLR2"] - sub["sc_Macro_FOLR2"].mean()) / sub["sc_Macro_FOLR2"].std()
        new_label = pd.Series("Macro_general", index=sub.index)
        # only-SPP1
        m_only_spp1  = spp1_pass & ~folr2_pass
        # only-FOLR2
        m_only_folr2 = ~spp1_pass & folr2_pass
        # both: pick larger z
        m_both = spp1_pass & folr2_pass
        new_label.loc[m_only_spp1]  = "Macro_SPP1"
        new_label.loc[m_only_folr2] = "Macro_FOLR2"
        if m_both.sum():
            new_label.loc[m_both] = np.where(spp1_z[m_both] >= folr2_z[m_both],
                                             "Macro_SPP1", "Macro_FOLR2")
        raw.obs.loc[is_general, "cell_type_fine"] = new_label.values
        print(f"   derived Macro_SPP1: {(new_label=='Macro_SPP1').sum()}, "
              f"Macro_FOLR2: {(new_label=='Macro_FOLR2').sum()}, "
              f"remaining Macro_general: {(new_label=='Macro_general').sum()}")

    # 7. Pull neu_subtype from luad_neutrophil_own_annotated.h5ad
    print("[7] merging neu_subtype ...")
    neu = load_obs_only(PROC / "luad_neutrophil_own_annotated.h5ad")
    neu = neu[["neu_subtype"]]
    # drop unclassified
    neu_keep = neu[neu["neu_subtype"].astype(str) != "Neu_unclassified"]
    common = raw.obs_names.intersection(neu_keep.index)
    print(f"   neutrophil cells overlap: {len(common)} (after drop unclassified)")
    raw.obs.loc[common, "cell_type_fine"] = neu_keep.loc[common, "neu_subtype"].astype(str).values

    # Also remove neutrophils without semantic subtype (but in myeloid as 'Neutrophil')
    is_neu_legacy = (raw.obs["cell_type_fine"].values == "Neutrophil")
    n_neu_legacy = int(is_neu_legacy.sum())
    if n_neu_legacy:
        print(f"   dropping {n_neu_legacy} Neutrophil cells without fine subtype")
        keep_idx = ~is_neu_legacy
        raw = raw[keep_idx].copy()
        gc.collect()

    # 8. Override Malignant epithelial cells
    print("[8] merging malignant flag ...")
    mal = load_obs_only(PROC / "luad_malignant_scored.h5ad")
    mal = mal[["malignant", "dominant_MP"]]
    mal_pos = mal[mal["malignant"].astype(str).str.lower().isin(["true","1","malignant"])]
    common = raw.obs_names.intersection(mal_pos.index)
    print(f"   malignant cells overlap: {len(common)}")
    raw.obs.loc[common, "cell_type_fine"] = "Malignant"

    # 9. Final cleanup: drop tiny categories, drop NaN/empty
    raw = raw[~raw.obs["cell_type_fine"].isna()].copy()
    raw = raw[raw.obs["cell_type_fine"].astype(str) != ""].copy()
    counts = raw.obs["cell_type_fine"].value_counts()
    drop_types = counts[counts < MIN_PER_TYPE].index.tolist()
    if drop_types:
        print(f"   dropping rare types (n<{MIN_PER_TYPE}): {drop_types}")
        raw = raw[~raw.obs["cell_type_fine"].isin(drop_types)].copy()

    # 10. Cap per cell type
    print(f"[10] capping each cell type at {MAX_PER_TYPE} cells ...")
    keep_idx_arr = []
    obs_index_arr = np.arange(raw.n_obs)
    for ct, sub_idx in raw.obs.groupby("cell_type_fine", observed=True).indices.items():
        if len(sub_idx) > MAX_PER_TYPE:
            choice = RNG.choice(sub_idx, MAX_PER_TYPE, replace=False)
        else:
            choice = sub_idx
        keep_idx_arr.append(choice)
    keep_idx_arr = np.sort(np.concatenate(keep_idx_arr))
    raw = raw[keep_idx_arr].copy()
    gc.collect()

    # Final counts report
    final_counts = raw.obs["cell_type_fine"].value_counts().sort_values(ascending=False)
    print("\n=== final cell_type_fine counts (after cap) ===")
    for ct, n in final_counts.items():
        print(f"   {ct}: {n}")
    print(f"\nTOTAL cells: {raw.n_obs}, genes: {raw.n_vars}")

    final_counts.to_frame("n_cells").to_csv(OUT / "cell_counts_report.csv")
    out_path = OUT / "unified_reference.h5ad"
    print(f"\n[save] {out_path}")
    raw.write_h5ad(str(out_path), compression="gzip")
    print(f"[save] file size: {out_path.stat().st_size/1e9:.2f} GB")

    # Sanity: verify integer counts preserved
    if sp.issparse(raw.X):
        sample = raw.X[: min(2000, raw.n_obs)].toarray().ravel()
    else:
        sample = np.asarray(raw.X[: min(2000, raw.n_obs)]).ravel()
    nz = sample[sample != 0]
    print(f"[sanity] post-save nonzero values: min={nz.min()}, max={nz.max()}, "
          f"all integer={bool(np.all(np.equal(np.mod(nz,1),0)))}")


if __name__ == "__main__":
    main()
