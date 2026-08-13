"""F2 — Batch of sensitivity analyses (Revision analysis:, 10, 13, 20 and 23).

S1 Sensitivity of spatial ROI z-thresholds (comments 13 and 23)
   Main text defines the NF-κB–neutrophil niche and tumour-intrinsic ROI at z>0.5; called this arbitrary.
   Recompute at z>0 / 0.25 / 0.5 / 0.75 / 1.0 and check whether key conclusions change with the threshold.

S2 Sensitivity of diffusion-pseudotime root cell
   Main text uses the highest-MP4 cell as root; this is biologically motivated but arbitrary.
   Resample roots from the top 10/50/100 MP4 cells and from the highest MP1/MP2 cell; check whether pseudotime ranks remain stable.

S3 Whether MP2 and MP3 co-occur in the same brain-metastasis cells
   Main text notes co-enrichment of both in brain metastases without stating same-cell co-occurrence.

S4 Direct association of Macro_SPP1 with MP3
   Fig. 4H lacks a direct Macro_SPP1–MP3 correlation, yet that axis is described as pro-EMT.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, os, warnings
from scipy import stats
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/sensitivity"
os.makedirs(OUT, exist_ok=True)

print("=== S1 Sensitivity of spatial ROI z-thresholds ===", flush=True)
a = ad.read_h5ad("${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad")
c = ad.read_h5ad("${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad")
com = a.obs_names.intersection(c.obs_names); a = a[com]; c = c[com]
z = lambda v: (np.asarray(v, float) - np.asarray(v, float).mean()) / (np.asarray(v, float).std(ddof=0) + 1e-12)
zn = z(a.obs["progeny_NFkB"].values) if "progeny_NFkB" in a.obs else None
neu_cols = [k for k in c.obsm["q05_cell_abundance_w_sf"].columns if "Neu_" in k]
zneu = z(c.obsm["q05_cell_abundance_w_sf"][neu_cols].sum(axis=1).values)
zmp3 = z(a.obs["MP3_score"].values)
zmal = z(c.obsm["q05_cell_abundance_w_sf"]["q05cell_abundance_w_sf_Malignant"].values)
X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
gi = {g: i for i, g in enumerate(a.var_names)}
AP1 = [g for g in ["JUNB", "FOS", "FOSB", "JUN", "ATF3", "NFKBIA", "IL1B"] if g in gi]

print(f"{'threshold':>6}{'stromal-immune ROI spots':>18}{'mean AP-1 Δ':>12}{'tumour-intrinsic ROI':>13}", flush=True)
rows = []
for th in [0.0, 0.25, 0.5, 0.75, 1.0]:
    r1 = (zn > th) & (zneu > th) if zn is not None else None
    r2 = (zmp3 > th) & (zmal > th)
    d1 = np.nan
    if r1 is not None and r1.sum() > 20:
        d1 = np.mean([X[r1, gi[g]].mean() - X[~r1, gi[g]].mean() for g in AP1])
    print(f"  {th:>6.2f}{int(r1.sum()) if r1 is not None else 0:>18}{d1:>12.3f}{int(r2.sum()):>13}", flush=True)
    rows.append(dict(threshold=th, n_niche=int(r1.sum()) if r1 is not None else 0,
                     ap1_delta=d1, n_tumor_roi=int(r2.sum())))
pd.DataFrame(rows).to_csv(f"{OUT}/S1_spatial_threshold.csv", index=False)

print("\n=== S3 Whether MP2 and MP3 co-occur in the same brain-metastasis cells ===", flush=True)
mal = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad")
o = mal.obs
brain = o.tissue_type.astype(str) == "Brain_Metastasis"
print(f"brain-metastasis malignant cells: {int(brain.sum())}", flush=True)
if brain.sum() > 50:
    b = o.loc[brain]
    print(f"dominant_MP distribution: {b.dominant_MP.value_counts().to_dict()}", flush=True)
    r, p = stats.spearmanr(b.MP2_score, b.MP3_score)
    ra, pa = stats.spearmanr(o.MP2_score, o.MP3_score)
    print(f"within brain mets, MP2 vs MP3 score Spearman ρ={r:+.3f} (p={p:.2g})", flush=True)
    print(f"all malignant cells ρ={ra:+.3f} (p={pa:.2g})", flush=True)
    hi2 = b.MP2_score > b.MP2_score.quantile(.75)
    hi3 = b.MP3_score > b.MP3_score.quantile(.75)
    obs_ = (hi2 & hi3).mean(); exp_ = hi2.mean() * hi3.mean()
    print(f"cells in the upper quartile of both: {obs_*100:.1f}% (independence expectation {exp_*100:.1f}%),"
          f"ratio {obs_/exp_:.2f}", flush=True)
    pd.DataFrame([dict(rho_brain=r, p_brain=p, rho_all=ra, obs_co=obs_, exp_co=exp_)]).to_csv(
        f"{OUT}/S3_mp2_mp3_cooccur.csv", index=False)

print("\n=== S4 Patient-level correlation of Macro_SPP1 with each MP ===", flush=True)
mye = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad")
if "myeloid_subtype" in mye.obs:
    mo = mye.obs
    key = "dataset" if "dataset" in mo else None
    mo = mo.assign(pt=mo["dataset"].astype(str) + "__" + mo["patient_id"].astype(str))
    o2 = o.assign(pt=o["dataset"].astype(str) + "__" + o["patient_id"].astype(str))
    sub_col = [c for c in ["dominant_macro", "myeloid_subtype"] if c in mo.columns][0]
    frac = (mo.groupby(["pt", sub_col]).size().unstack(fill_value=0)
            .pipe(lambda d: d.div(d.sum(1), axis=0)))
    mpm = o2.groupby("pt")[["MP1_score", "MP2_score", "MP3_score", "MP4_score"]].mean()
    j = frac.join(mpm, how="inner")
    tgt = [c for c in frac.columns if "SPP1" in str(c)]
    print(f"paired patients: {len(j)}; SPP1-related columns: {tgt}", flush=True)
    for t in tgt:
        line = f"  {t}: "
        for mp in ["MP1_score", "MP2_score", "MP3_score", "MP4_score"]:
            r, p = stats.spearmanr(j[t], j[mp])
            line += f"{mp[:3]} ρ={r:+.3f}({'*' if p<0.05 else ' '}p={p:.2g})  "
        print(line, flush=True)
    j.to_csv(f"{OUT}/S4_macro_mp_patient.csv")
print(f"\nResults written to {OUT}/", flush=True)
