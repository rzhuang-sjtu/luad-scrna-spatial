"""Step 14: re-export Fig 2 panels at consensus-cNMF program level (15 programs).

Cuts the existing step6 770-GEP dendrogram at 15 clusters, derives consensus
program signatures by mean GEP gene weight, then exports:

  1. cnmf_consensus_corr.csv         15×15 Spearman over consensus signatures
  2. cnmf_consensus_mp_annotation.csv  cluster → MP (majority vote + purity)
  3. monocle_trajectory.csv.gz       cells × DC1/DC2/pseudotime/MP
  4. gavish_cosine_similarity_cnmf.csv  15 × 41 cosine similarity vs Gavish

Inputs:
  - ~/luad/data/cnmf_output/gep_pool_tpm.csv     9772 genes × 770 GEPs
  - ~/luad/data/cnmf_output/gep_top100_genes.csv 770 × top100 gene
  - ~/luad/data/cnmf_output/gep_metadata.csv
  - ~/luad/results/step6_linkage.npz             hierarchical linkage
  - ~/luad/results/step6_gep_mp_assignment.csv   per-GEP MP label
  - ~/luad/data/processed/luad_malignant_scored.h5ad
  - ~/luad/results/step10b_pseudotime.csv.gz     existing DPT pseudotime
  - ~/luad/data/reference/gavish2023_MPs.csv     41 Gavish MPs × 50 genes
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

DATA = Path.home() / "luad/data"
RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig2")
N_CONSENSUS = 15
TOPN_GENES = 50  # for Gavish cosine


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()

    # 1. Load 770-GEP pool + linkage + per-GEP MP assignment
    log("loading GEP pool + step6 linkage")
    gep_pool = pd.read_csv(DATA/"cnmf_output/gep_pool_tpm.csv", index_col=0)
    log(f"  GEP pool TPM: {gep_pool.shape}  (genes × GEPs)")

    gep_meta = pd.read_csv(DATA/"cnmf_output/gep_metadata.csv")
    log(f"  GEP metadata: {gep_meta.shape}")

    mp_assign = pd.read_csv(RES/"step6_gep_mp_assignment.csv")
    log(f"  GEP→MP: {mp_assign.shape}")

    Z = np.load(RES/"step6_linkage.npz")
    log(f"  linkage keys: {list(Z.keys())}")
    linkage = Z["linkage"] if "linkage" in Z.files else Z[list(Z.keys())[0]]

    from scipy.cluster.hierarchy import fcluster
    cluster_labels = fcluster(linkage, t=N_CONSENSUS, criterion="maxclust")
    log(f"  cut to {len(np.unique(cluster_labels))} clusters")

    # leaves order in linkage matches columns of gep_pool? Need to verify.
    # The dendrogram order csv stores it. For fcluster, indices follow the
    # original input order (here, the GEP columns of gep_pool which match
    # gep_meta).
    if len(cluster_labels) != gep_pool.shape[1]:
        log(f"  WARN: cluster_labels {len(cluster_labels)} != GEPs {gep_pool.shape[1]}")
    gep_meta = gep_meta.copy()
    gep_meta["consensus_cluster"] = cluster_labels

    # Cluster ordering: rename clusters by descending size, prefix cNMF_
    sizes = gep_meta["consensus_cluster"].value_counts()
    rename = {old: f"cNMF_{i+1}" for i, old in enumerate(sizes.index)}
    gep_meta["program_id"] = gep_meta["consensus_cluster"].map(rename)
    log("  consensus cluster sizes:")
    log(gep_meta["program_id"].value_counts().sort_index().to_string())

    # 2. Consensus signatures = mean GEP TPM within each cluster
    log("computing consensus signatures (mean GEP TPM per cluster)")
    gep_to_pid = dict(zip(gep_meta["gep_id"], gep_meta["program_id"]))
    pids_in_order = sorted(set(gep_to_pid.values()),
                            key=lambda x: int(x.split("_")[1]))
    consensus_sig = pd.DataFrame(index=gep_pool.index, columns=pids_in_order,
                                  dtype=np.float32)
    for pid in pids_in_order:
        members = [g for g, p in gep_to_pid.items() if p == pid]
        consensus_sig[pid] = gep_pool[members].mean(axis=1).astype(np.float32)
    log(f"  consensus_sig shape: {consensus_sig.shape}")

    # OUTPUT 1: consensus correlation (Spearman)
    log("Spearman correlation across consensus signatures")
    # Use rank-then-pearson for speed (==Spearman)
    ranks = consensus_sig.rank(axis=0)
    corr = ranks.corr(method="pearson")  # rank-pearson = spearman
    corr.index.name = "program_id"
    corr.to_csv(FIG/"cnmf_consensus_corr.csv")
    log(f"  cnmf_consensus_corr.csv written {corr.shape}")

    # OUTPUT 2: cluster → MP majority vote + purity
    log("consensus → MP annotation")
    mp_lookup = dict(zip(mp_assign["gep_id"], mp_assign["MP"])) \
                if "MP" in mp_assign.columns else None
    if mp_lookup is None:
        # try alternate column name
        log(f"  mp_assign cols: {list(mp_assign.columns)}")
        # Fall back: pick the second column as MP label
        mp_col = [c for c in mp_assign.columns if c.lower() != "gep_id"][0]
        mp_lookup = dict(zip(mp_assign["gep_id"], mp_assign[mp_col]))
    gep_meta["MP"] = gep_meta["gep_id"].map(mp_lookup).fillna("Unassigned")

    rows = []
    for pid in pids_in_order:
        sub = gep_meta[gep_meta["program_id"] == pid]
        vc = sub["MP"].value_counts()
        top_mp = vc.index[0]
        purity = vc.iloc[0] / vc.sum()
        rows.append({"program_id": pid, "n_GEPs": len(sub),
                     "MP": top_mp, "purity": float(purity),
                     "MP_breakdown": dict(vc.to_dict())})
    ann = pd.DataFrame(rows)
    ann_out = ann[["program_id", "MP", "n_GEPs", "purity"]].copy()
    ann_out["MP_breakdown"] = ann["MP_breakdown"].apply(
        lambda d: ";".join(f"{k}:{v}" for k, v in d.items()))
    ann_out.to_csv(FIG/"cnmf_consensus_mp_annotation.csv", index=False)
    log(f"  cnmf_consensus_mp_annotation.csv:")
    log(ann_out.to_string(index=False))

    # OUTPUT 4: Gavish cosine similarity (consensus × Gavish)
    log("computing Gavish cosine similarity (binary signature vectors)")
    gavish = pd.read_csv(DATA/"reference/gavish2023_MPs.csv")
    gav_mps = {mp: set(df.gene) for mp, df in gavish.groupby("MP", sort=False)}

    # Top-50 genes per consensus program (by mean TPM in consensus_sig)
    consensus_top = {}
    for pid in pids_in_order:
        top = consensus_sig[pid].sort_values(ascending=False).head(TOPN_GENES)
        consensus_top[pid] = set(top.index.tolist())

    # Universe = union of all genes appearing in any signature (both sides)
    universe = set()
    for s in consensus_top.values(): universe |= s
    for s in gav_mps.values(): universe |= s
    universe = sorted(universe)
    log(f"  universe size: {len(universe)}")
    universe_idx = {g: i for i, g in enumerate(universe)}

    def to_vec(gene_set):
        v = np.zeros(len(universe), dtype=np.float32)
        for g in gene_set:
            i = universe_idx.get(g)
            if i is not None: v[i] = 1.0
        return v

    consensus_vecs = {pid: to_vec(s) for pid, s in consensus_top.items()}
    gav_vecs = {mp: to_vec(s) for mp, s in gav_mps.items()}

    # Cosine similarity
    cos_mat = np.zeros((len(pids_in_order), len(gav_vecs)), dtype=np.float32)
    for i, pid in enumerate(pids_in_order):
        a = consensus_vecs[pid]; na = np.linalg.norm(a)
        for j, gn in enumerate(gav_vecs):
            b = gav_vecs[gn]; nb = np.linalg.norm(b)
            cos_mat[i, j] = float(np.dot(a, b) / (na * nb)) if (na*nb) > 0 else 0.0
    cos_df = pd.DataFrame(cos_mat, index=pids_in_order,
                            columns=list(gav_vecs.keys()))
    cos_df.index.name = "program_id"
    cos_df.to_csv(FIG/"gavish_cosine_similarity_cnmf.csv")
    log(f"  gavish_cosine_similarity_cnmf.csv written {cos_df.shape}")
    # quick sanity: top match per consensus
    log("  top Gavish match per consensus program (by cosine):")
    for pid in pids_in_order:
        top = cos_df.loc[pid].sort_values(ascending=False).head(1)
        log(f"    {pid:8s} -> {top.index[0][:35]:35s} cos={top.iloc[0]:.3f}")

    # OUTPUT 3: Monocle-style trajectory CSV (DC1/DC2/pseudotime)
    log("loading malignant h5ad for diffmap")
    a = sc.read_h5ad(DATA/"processed/luad_malignant_scored.h5ad")
    log(f"  malignant: {a.shape}")

    if "X_diffmap" not in a.obsm:
        log("  rebuilding neighbors on X_pca_mal_harmony + diffmap")
        sc.pp.neighbors(a, use_rep="X_pca_mal_harmony", n_neighbors=30,
                         random_state=0)
        sc.tl.diffmap(a, n_comps=5, random_state=0)
    diff = a.obsm["X_diffmap"]
    # First diffusion component is trivial; use DC1/DC2 (cols 1, 2)
    dc1 = diff[:, 1]
    dc2 = diff[:, 2]
    log(f"  diffmap shape: {diff.shape}; using cols 1,2 as DC1/DC2")

    # Pseudotime from step10b
    pt_df = pd.read_csv(RES/"step10b_pseudotime.csv.gz")[["barcode","pseudotime"]]
    pt_lookup = dict(zip(pt_df["barcode"], pt_df["pseudotime"]))

    out = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "Component_1": dc1, "Component_2": dc2,
        "pseudotime": [pt_lookup.get(bc, np.nan) for bc in a.obs.index],
        "dominant_MP": a.obs["dominant_MP"].astype(str).values,
    })
    out.to_csv(FIG/"monocle_trajectory.csv.gz", index=False, compression="gzip")
    log(f"  monocle_trajectory.csv.gz written {out.shape}")

    # Summary
    log(f"\nALL 4 outputs written to {FIG}:")
    log("  - cnmf_consensus_corr.csv")
    log("  - cnmf_consensus_mp_annotation.csv")
    log("  - gavish_cosine_similarity_cnmf.csv")
    log("  - monocle_trajectory.csv.gz")
    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
