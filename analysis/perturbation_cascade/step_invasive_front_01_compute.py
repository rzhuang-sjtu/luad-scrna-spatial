"""Three-ROI analysis: 3rd ROI definition + co-occurrence + cross-ROI gene tests.

The third region is labelled Compositional-mixing in the figures; the column
names below keep the older invasive_front spelling, which is what the stored
result tables contain.

Outputs to ${DATA_ROOT}/ST/results/step_invasive_front/:
  - cooccurrence_long.csv   :  per-section P(target | anchor at distance d)
  - roi_overlap.csv         :  3×3 overlap matrix per cohort
  - gene_by_roi_stats.csv   :  7 genes × 3 ROI × 2 cohorts MW test
  - spots_{cohort}.csv      :  spot-level x,y,z_mal,z_neu,z_mp3,3 ROI flags
  - representative_sections.csv : section with largest invasive-front ROI per cohort

ROI definitions:
  - stromal_immune    = original ROI from cohort_with_roi.h5ad ("roi" col)
                          ≈ z(NFkB)>0.5 AND z(Neu_total)>0.5
  - tumor_intrinsic   = z(Malignant)>0.5 AND z(MP3_score)>0.5  (per-cohort z)
  - invasive_front    = z(Malignant)∈(0,1) AND z(Neu_total)∈(0,1)
                          (both elevated but neither extreme — transition zone)

Co-occurrence: per-spot 4-class assignment (argmax of z over Mal/Neu/SPP1/Fib;
"Other" if max z < 0.5). For each pair (Mal-Neu, Mal-SPP1, Mal-Fib), compute
P(target | anchor at distance d) via KDTree. Distance bins: 0-2500 px step 100.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KDTree
from scipy import stats
warnings.filterwarnings("ignore")

OUT = Path("${DATA_ROOT}/ST/results/step_invasive_front")
OUT.mkdir(parents=True, exist_ok=True)

GENES = ["SEC61G", "SRSF9", "ANGPTL4", "OSM", "IL1B", "ATF3", "FOSB"]
NEU_SUBTYPES = ["Neu_Inflammatory", "Neu_Angiogenic", "Neu_Metastatic",
                "Neu_ECM_remodeling", "Neu_OSM_priming", "Neu_OSM_low",
                "Neu_IFN_response"]
COOCCUR_PAIRS = [("Malignant", "Neutrophil"),
                 ("Malignant", "Macro_SPP1"),
                 ("Malignant", "Fibroblast")]
INTERVAL_MAX  = 2500.0   # full-res pixel distance, ~770 μm at standard Visium scaling
INTERVAL_STEP = 100.0

DATASETS = {
    "E-MTAB-13530": dict(
        c2l="${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad",
        roi="${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad"),
    "Okamura":      dict(
        c2l="${DATA_ROOT}/ST/results/step09_okamura_validation/all_sections_c2l.h5ad",
        roi="${DATA_ROOT}/ST/results/step09_okamura_validation/cohort_with_roi.h5ad"),
}

def z(x):
    x = np.asarray(x, dtype=float)
    s = x.std(ddof=0)
    return (x - x.mean()) / s if s > 0 else np.zeros_like(x)

def cooccur_section(coords, cls, intervals, anchor, target):
    tree = KDTree(coords)
    n_bins = len(intervals) - 1
    counts_t = np.zeros(n_bins)
    counts_a = np.zeros(n_bins)
    is_anchor = (cls == anchor)
    is_target = (cls == target)
    if not is_anchor.any():
        return np.full(n_bins, np.nan), np.zeros(n_bins, dtype=int)
    inds, dists = tree.query_radius(coords[is_anchor],
                                    r=intervals[-1], return_distance=True)
    anchor_idxs = np.where(is_anchor)[0]
    for k, (j_arr, d_arr) in enumerate(zip(inds, dists)):
        i = anchor_idxs[k]
        m = (j_arr != i) & (d_arr > 0)
        if not m.any(): continue
        j_arr, d_arr = j_arr[m], d_arr[m]
        bin_idx = np.searchsorted(intervals, d_arr, side="right") - 1
        valid = (bin_idx >= 0) & (bin_idx < n_bins)
        bin_idx = bin_idx[valid]
        np.add.at(counts_a, bin_idx, 1)
        np.add.at(counts_t, bin_idx[is_target[j_arr[valid]]], 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        prob = np.where(counts_a > 0, counts_t / counts_a, np.nan)
    return prob, counts_a.astype(int)

all_cooccur, all_overlap, all_gene_stats, all_repsec = [], [], [], []

for label, cfg in DATASETS.items():
    print(f"\n=== {label} ===")
    ad_c2l = sc.read_h5ad(cfg["c2l"])
    ad_roi = sc.read_h5ad(cfg["roi"])
    common = ad_c2l.obs_names.intersection(ad_roi.obs_names)
    ad_c2l = ad_c2l[common].copy()
    ad_roi = ad_roi[common].copy()
    print(f"  spots: {len(common)}; sections: {ad_c2l.obs['sample'].nunique()}")

    abund = ad_c2l.obsm["q05_cell_abundance_w_sf"]
    mal  = abund["q05cell_abundance_w_sf_Malignant"].values
    sppl = abund["q05cell_abundance_w_sf_Macro_SPP1"].values
    fib  = abund["q05cell_abundance_w_sf_Fibroblast"].values
    neu_cols = [f"q05cell_abundance_w_sf_{s}" for s in NEU_SUBTYPES]
    neu_total = abund[neu_cols].sum(axis=1).values
    mp3 = ad_roi.obs["MP3_score"].values.astype(float)

    z_mal, z_neu, z_sppl, z_fib, z_mp3 = z(mal), z(neu_total), z(sppl), z(fib), z(mp3)

    # 4-class spot assignment
    Z = np.column_stack([z_mal, z_neu, z_sppl, z_fib])
    cls_idx = np.argmax(Z, axis=1)
    cls_names = np.array(["Malignant", "Neutrophil", "Macro_SPP1", "Fibroblast"])
    cls = cls_names[cls_idx]
    cls[Z.max(axis=1) < 0.5] = "Other"
    counts = pd.Series(cls).value_counts()
    print(f"  4-class assignment: {counts.to_dict()}")

    # ROIs (per-cohort z, applied globally)
    orig_roi   = ad_roi.obs["roi"].values.astype(bool)
    tumor_intr = (z_mal > 0.5) & (z_mp3 > 0.5)
    invasive   = (z_mal > 0) & (z_neu > 0) & (z_mal < 1) & (z_neu < 1)
    print(f"  ROI sizes: stromal-immune={orig_roi.sum()}, tumor-intrinsic={tumor_intr.sum()}, invasive-front={invasive.sum()}")

    rois = {"stromal_immune": orig_roi, "tumor_intrinsic": tumor_intr, "invasive_front": invasive}
    for a in rois:
        for b in rois:
            ov = int((rois[a] & rois[b]).sum())
            all_overlap.append({"dataset": label, "ROI_A": a, "ROI_B": b,
                                "n_overlap": ov,
                                "n_A": int(rois[a].sum()),
                                "n_B": int(rois[b].sum()),
                                "jaccard": ov / max(int((rois[a] | rois[b]).sum()), 1)})

    # 7 genes × 3 ROIs MW test
    for g in GENES:
        if g not in ad_roi.var_names:
            print(f"  WARN gene not in data: {g}")
            continue
        e = ad_roi[:, g].X
        e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
        for roi_name, mask in rois.items():
            a_e, b_e = e[mask], e[~mask]
            if len(a_e) < 5 or len(b_e) < 5: continue
            u, p = stats.mannwhitneyu(a_e, b_e, alternative="two-sided")
            all_gene_stats.append({
                "dataset": label, "gene": g, "roi_type": roi_name,
                "n_roi": int(mask.sum()), "n_non": int((~mask).sum()),
                "mean_roi": float(a_e.mean()), "mean_non": float(b_e.mean()),
                "delta_roi_minus_non": float(a_e.mean() - b_e.mean()),
                "mw_p": float(p),
            })

    # Per-section co-occurrence
    intervals = np.arange(0, INTERVAL_MAX + INTERVAL_STEP, INTERVAL_STEP)
    for s in sorted(ad_c2l.obs["sample"].unique()):
        m = ad_c2l.obs["sample"].values == s
        coords = ad_c2l.obsm["spatial"][m]
        cls_s = cls[m]
        for anchor, target in COOCCUR_PAIRS:
            prob, n_a = cooccur_section(coords, cls_s, intervals, anchor, target)
            for i in range(len(prob)):
                all_cooccur.append({
                    "dataset": label, "sample": s,
                    "anchor": anchor, "target": target,
                    "dist_low": float(intervals[i]),
                    "dist_high": float(intervals[i+1]),
                    "dist_mid": float((intervals[i] + intervals[i+1]) / 2),
                    "prob": float(prob[i]) if not np.isnan(prob[i]) else np.nan,
                    "n_anchor_neighbors_in_bin": int(n_a[i]),
                })
    print(f"  cooccurrence done")

    # Spots CSV (for plotting H&E overlay)
    df_spots = pd.DataFrame({
        "obs_name": ad_c2l.obs_names,
        "sample":  ad_c2l.obs["sample"].values,
        "x":       ad_c2l.obsm["spatial"][:, 0],
        "y":       ad_c2l.obsm["spatial"][:, 1],
        "z_mal":   z_mal, "z_neu": z_neu, "z_mp3": z_mp3,
        "stromal_immune_roi":  orig_roi,
        "tumor_intrinsic_roi": tumor_intr,
        "invasive_front_roi":  invasive,
        "spot_class": cls,
        "dataset": label,
    })
    df_spots.to_csv(OUT / f"spots_{label}.csv", index=False)

    # Representative section: largest invasive_front count
    inv_per_sec = pd.Series(invasive).groupby(ad_c2l.obs["sample"].values).sum()
    rep = inv_per_sec.idxmax()
    all_repsec.append({"dataset": label, "section": str(rep),
                       "n_invasive_front": int(inv_per_sec.max()),
                       "n_stromal_immune_in_section": int(orig_roi[ad_c2l.obs["sample"].values == rep].sum()),
                       "n_tumor_intrinsic_in_section": int(tumor_intr[ad_c2l.obs["sample"].values == rep].sum()),
                       "n_total_spots_in_section": int((ad_c2l.obs["sample"].values == rep).sum())})
    print(f"  representative section: {rep}")

pd.DataFrame(all_cooccur).to_csv(OUT / "cooccurrence_long.csv", index=False)
pd.DataFrame(all_overlap).to_csv(OUT / "roi_overlap.csv", index=False)
pd.DataFrame(all_gene_stats).to_csv(OUT / "gene_by_roi_stats.csv", index=False)
pd.DataFrame(all_repsec).to_csv(OUT / "representative_sections.csv", index=False)

# pretty print summary
print("\n" + "=" * 80)
print(" 3-ROI overlap (jaccard)")
print("=" * 80)
ov_df = pd.DataFrame(all_overlap)
for ds in ov_df["dataset"].unique():
    sub = ov_df[ov_df["dataset"] == ds].pivot(index="ROI_A", columns="ROI_B", values="jaccard")
    print(f"\n{ds}:"); print(sub.round(3).to_string())

print("\n" + "=" * 80)
print(" 7 genes × 3 ROI × 2 cohorts (delta and -log10(p))")
print("=" * 80)
gs = pd.DataFrame(all_gene_stats)
gs["star"] = gs["mw_p"].apply(lambda p: ("****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"))
piv = gs.pivot_table(index=["gene"], columns=["dataset", "roi_type"],
                     values=["delta_roi_minus_non"]).round(3)
print(piv.to_string())
print(f"\nAll outputs in: {OUT}")
