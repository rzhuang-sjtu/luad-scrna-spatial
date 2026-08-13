"""F3 — Three small sensitivity analyses.

S2 Sensitivity of diffusion-pseudotime root cell (sensitivity to the pseudotime root cell)
   Main text uses the single highest-MP4 cell as root; this is biologically motivated but arbitrary.
   Resample roots at random from the top 10 / 50 / 100 MP4 cells, and from the highest MP1 / MP2 cell,
   and check stability of cell pseudotime ranks (Spearman across schemes).

S5 cNMF K_mp choice and MP5 exclusion criteria (comments 8, 17)
   Main text scans K_mp∈[3,8] by silhouette plus patient-mixing entropy, selects K_mp=5, then drops MP5.
   Rerun that scan and report silhouette and per-cluster patient entropy at each K to justify the choice.

S6 Sensitivity of restricting the Salcher neutrophil panel to LUAD
   The Salcher reference atlas is NSCLC (including squamous); its use is questionable for pure LUAD.
   Compare label transfer when training on all reference cells vs adenocarcinoma samples only.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, os, warnings
from scipy import stats
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/sensitivity"
os.makedirs(OUT, exist_ok=True)

print("=== S2 Sensitivity of diffusion-pseudotime root cell ===", flush=True)
mal = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad")
if "X_diffmap" in mal.obsm:
    DM = np.asarray(mal.obsm["X_diffmap"])[:, 1:5]      # first 4 non-trivial components
    print(f"using existing diffmap, {DM.shape}", flush=True)
else:
    import scanpy as sc
    key = [k for k in ["X_pca_mal_harmony", "X_pca_harmony", "X_pca"] if k in mal.obsm][0]
    print(f"no diffmap; compute from {key} ...", flush=True)
    sc.pp.neighbors(mal, use_rep=key, n_neighbors=30)
    sc.tl.diffmap(mal, n_comps=5)
    DM = np.asarray(mal.obsm["X_diffmap"])[:, 1:5]

o = mal.obs
rng = np.random.default_rng(20260805)

def pseudotime_from(root_idx):
    """Diffusion distance to the root cell as a simplified pseudotime proxy"""
    d = np.linalg.norm(DM - DM[root_idx], axis=1)
    return d

schemes = {}
schemes["main text: highest MP4"] = [int(np.argmax(o.MP4_score.values))]
for k in [10, 50, 100]:
    top = np.argsort(-o.MP4_score.values)[:k]
    schemes[f"random among top {k} MP4"] = list(rng.choice(top, 5, replace=False))
schemes["highest MP1"] = [int(np.argmax(o.MP1_score.values))]
schemes["highest MP2"] = [int(np.argmax(o.MP2_score.values))]

ref = pseudotime_from(schemes["main text: highest MP4"][0])
rows = []
for name, roots in schemes.items():
    rhos = [stats.spearmanr(ref, pseudotime_from(r))[0] for r in roots]
    # Whether mean-MP pseudotime order is preserved (main text: MP4 lowest < MP1 ≈ MP3 < MP2)
    pt = pseudotime_from(roots[0])
    order = o.assign(pt=pt).groupby("dominant_MP").pt.mean().sort_values()
    rows.append(dict(scheme=name, rho_vs_ref=np.mean(rhos),
                     mp_order="<".join(order.index.tolist())))
    print(f"{name:<18} Spearman vs main-text ranks={np.mean(rhos):+.3f}   MP order:"
          f"{' < '.join(order.index.tolist())}", flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/S2_root_sensitivity.csv", index=False)

print("\n=== S5 cNMF K_mp choice and MP5 exclusion ===", flush=True)
gp = "${PROJECT_ROOT}/data/cnmf_output/gep_pool_zscore.csv"
if os.path.exists(gp):
    Z = pd.read_csv(gp, index_col=0)
    print(f"GEP pool {Z.shape} (genes x GEPs)", flush=True)
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    from sklearn.metrics import silhouette_score
    C = Z.corr(method="spearman")
    D = 1 - C.values
    np.fill_diagonal(D, 0)
    L = linkage(squareform(D, checks=False), method="average")
    meta = pd.read_csv("${PROJECT_ROOT}/data/cnmf_output/gep_metadata.csv")
    pat = meta.set_index(meta.columns[0]).iloc[:, 0] if len(meta.columns) > 1 else None
    print(f"{'K_mp':>5}{'silhouette':>12}{'min cluster':>8}{'cluster sizes':>28}", flush=True)
    for k in range(3, 9):
        lab = fcluster(L, k, criterion="maxclust")
        s = silhouette_score(D, lab, metric="precomputed")
        sz = pd.Series(lab).value_counts().sort_index()
        print(f"  {k:>5}{s:>12.4f}{sz.min():>8}   {list(sz.values)}", flush=True)
else:
    print("GEP pool file not found; skip", flush=True)
print(f"\nResults written to {OUT}/", flush=True)
