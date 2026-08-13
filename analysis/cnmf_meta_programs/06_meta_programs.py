"""
Step 6: GEP pool -> Spearman clustering -> meta-programs (MPs).

Pipeline:
    1. Load data/cnmf_output/gep_pool_zscore.csv (9881 genes x 799 GEPs)
    2. Pairwise Spearman correlation (GEP x GEP)
    3. Hierarchical clustering (average linkage on 1 - Spearman)
    4. Explore K_mp = 3..8; report for each:
        - silhouette on precomputed distance
        - n_patients per MP (patient mixing)
        - min/median MP size
    5. Pick K_mp (default 5, overridable) and emit:
        - step6_gep_correlation.csv           (799x799)
        - step6_linkage.npz                   (scipy Z matrix + leaf order)
        - step6_mp_scan.csv                   (K_mp x metrics table)
        - step6_gep_mp_assignment.csv         (one row per GEP)
        - step6_mp_signatures_top100.csv      (ranked genes per MP)
        - step6_mp_patient_mixing.csv         (n_GEPs, n_patients, entropy, H_norm)
        - step6_mp_summary.md                 (human-readable summary)
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import rankdata
from sklearn.metrics import silhouette_score

ROOT = Path("${PROJECT_ROOT}")
CNMF_OUTPUT = ROOT / "data" / "cnmf_output"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

POOL_CSV = CNMF_OUTPUT / "gep_pool_zscore.csv"
META_CSV = CNMF_OUTPUT / "gep_metadata.csv"
TPM_CSV = CNMF_OUTPUT / "gep_pool_tpm.csv"
TOP100_CSV = CNMF_OUTPUT / "gep_top100_genes.csv"

OUT_CORR = RESULTS / "step6_gep_correlation.csv"
OUT_LINK = RESULTS / "step6_linkage.npz"
OUT_SCAN = RESULTS / "step6_mp_scan.csv"
OUT_ASSIGN = RESULTS / "step6_gep_mp_assignment.csv"
OUT_SIGS = RESULTS / "step6_mp_signatures_top100.csv"
OUT_MIX = RESULTS / "step6_mp_patient_mixing.csv"
OUT_SUMMARY = RESULTS / "step6_mp_summary.md"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def spearman_corr(pool: pd.DataFrame) -> pd.DataFrame:
    """Rank each column then Pearson = Spearman. pool: genes x GEPs."""
    log(f"Ranking columns of pool ({pool.shape[0]} genes x {pool.shape[1]} GEPs) ...")
    ranks = np.empty_like(pool.values, dtype=np.float32)
    for j in range(pool.shape[1]):
        ranks[:, j] = rankdata(pool.values[:, j]).astype(np.float32)
    log("Computing Pearson on ranks ...")
    corr = np.corrcoef(ranks, rowvar=False).astype(np.float32)
    return pd.DataFrame(corr, index=pool.columns, columns=pool.columns)


def patient_mixing(assign: pd.DataFrame) -> pd.DataFrame:
    """Per-MP: n_GEPs, n_patients, Shannon entropy over patient GEP counts, H/log(n_patients)."""
    total_patients = assign["patient_key"].nunique()
    out = []
    for mp, sub in assign.groupby("MP"):
        n_geps = len(sub)
        pats = sub["patient_key"].value_counts()
        p = (pats / pats.sum()).values
        H = float(-(p * np.log(p + 1e-12)).sum())
        H_norm = H / np.log(total_patients) if total_patients > 1 else 0.0
        out.append({
            "MP": mp,
            "n_GEPs": n_geps,
            "n_patients": int(len(pats)),
            "entropy": H,
            "H_norm": H_norm,
        })
    return pd.DataFrame(out).sort_values("MP").reset_index(drop=True)


def mp_signatures(pool: pd.DataFrame, assign: pd.DataFrame, n_top: int = 100) -> pd.DataFrame:
    """Consensus signature per MP: average GEP z-score across genes, rank, take top.

    pool: genes x GEPs   (gep_pool_zscore)
    assign: gep_id, MP
    """
    rows = []
    for mp, sub in assign.groupby("MP"):
        member_geps = sub["gep_id"].tolist()
        mean_scores = pool[member_geps].mean(axis=1)
        top = mean_scores.sort_values(ascending=False).head(n_top)
        for rank, (gene, score) in enumerate(top.items(), start=1):
            rows.append({
                "MP": mp,
                "rank": rank,
                "gene": gene,
                "consensus_score": float(score),
                "n_GEPs_in_MP": len(member_geps),
            })
    return pd.DataFrame(rows)


def scan_k(D_cond: np.ndarray, D_square: np.ndarray, Z: np.ndarray,
           meta: pd.DataFrame, k_values) -> tuple[pd.DataFrame, dict]:
    rows = []
    per_k_assign = {}
    for k in k_values:
        labels = fcluster(Z, t=k, criterion="maxclust")
        sizes = pd.Series(labels).value_counts().sort_index()
        # silhouette on precomputed distance
        if len(set(labels)) > 1:
            sil = float(silhouette_score(D_square, labels, metric="precomputed"))
        else:
            sil = np.nan
        assign = meta.copy()
        assign["MP"] = labels
        mix = patient_mixing(assign)
        rows.append({
            "k": k,
            "silhouette": sil,
            "min_size": int(sizes.min()),
            "max_size": int(sizes.max()),
            "median_size": float(sizes.median()),
            "min_n_patients": int(mix["n_patients"].min()),
            "median_n_patients": float(mix["n_patients"].median()),
            "min_H_norm": float(mix["H_norm"].min()),
            "mean_H_norm": float(mix["H_norm"].mean()),
        })
        per_k_assign[k] = assign
    return pd.DataFrame(rows), per_k_assign


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=None,
                    help="Fixed K for final MP assignment. "
                         "If omitted, pick max silhouette in [4,6].")
    ap.add_argument("--linkage", default="average",
                    choices=["average", "complete", "ward"])
    ap.add_argument("--scan-k", default="3,4,5,6,7,8")
    args = ap.parse_args()

    log(f"Reading {POOL_CSV} ...")
    pool = pd.read_csv(POOL_CSV, index_col=0)
    log(f"Pool shape: {pool.shape}")
    meta = pd.read_csv(META_CSV)
    log(f"Metadata: {len(meta)} GEPs, {meta['patient_key'].nunique()} patients")

    corr = spearman_corr(pool)
    corr.to_csv(OUT_CORR, float_format="%.4f")
    log(f"Wrote {OUT_CORR} ({OUT_CORR.stat().st_size/1024/1024:.1f} MB)")

    # Distance + linkage
    D = (1.0 - corr.values).astype(np.float64)
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    D = np.clip(D, 0.0, 2.0)
    D_cond = squareform(D, checks=False)
    log(f"Linkage (method={args.linkage}) ...")
    if args.linkage == "ward":
        # Ward expects squared Euclidean-style. Feed condensed distance.
        Z = linkage(D_cond, method="ward")
    else:
        Z = linkage(D_cond, method=args.linkage)

    leaf_order = leaves_list(Z)
    gep_ids = corr.index.to_numpy()
    np.savez_compressed(OUT_LINK, Z=Z, leaf_order=leaf_order, gep_ids=gep_ids)
    log(f"Wrote {OUT_LINK}")

    # Scan K
    k_list = [int(x) for x in args.scan_k.split(",")]
    log(f"Scanning K_mp = {k_list} ...")
    scan_df, per_k_assign = scan_k(D_cond, D, Z, meta, k_list)
    scan_df.to_csv(OUT_SCAN, index=False)
    log("Scan results:")
    log(scan_df.to_string(index=False))
    log(f"Wrote {OUT_SCAN}")

    # Pick K
    if args.k is not None:
        k_final = args.k
        reason = f"fixed by --k={args.k}"
    else:
        # Prefer best silhouette in [4,6]; fallback full range.
        core = scan_df[scan_df["k"].isin([4, 5, 6])]
        if len(core) > 0 and core["silhouette"].notna().any():
            k_final = int(core.loc[core["silhouette"].idxmax(), "k"])
            reason = f"max silhouette in K=[4,6]"
        else:
            k_final = int(scan_df.loc[scan_df["silhouette"].idxmax(), "k"])
            reason = "max silhouette overall"
    log(f"Selected K_mp = {k_final}  ({reason})")

    assign = per_k_assign[k_final].copy()
    mp_counts = assign["MP"].value_counts().sort_index()
    # Relabel MPs by cluster size desc (MP1 = largest)
    size_order = mp_counts.sort_values(ascending=False).index.tolist()
    rename = {old: f"MP{new}" for new, old in enumerate(size_order, start=1)}
    assign["MP"] = assign["MP"].map(rename)
    assign = assign[["gep_id", "patient_key", "dataset", "k_total",
                     "program_idx", "n_cells", "stability", "MP"]]
    assign.to_csv(OUT_ASSIGN, index=False)
    log(f"Wrote {OUT_ASSIGN}")

    mix = patient_mixing(assign)
    mix.to_csv(OUT_MIX, index=False)
    log("Patient mixing:")
    log(mix.to_string(index=False))

    sigs = mp_signatures(pool, assign, n_top=100)
    sigs.to_csv(OUT_SIGS, index=False)
    log(f"Wrote {OUT_SIGS}")

    with open(OUT_SUMMARY, "w") as f:
        f.write(f"# Step 6 meta-program summary\n\n")
        f.write(f"- Pool: {pool.shape[0]} genes x {pool.shape[1]} GEPs\n")
        f.write(f"- Linkage: {args.linkage}\n")
        f.write(f"- Selected K_mp = **{k_final}** ({reason})\n\n")
        f.write("## Scan of K_mp\n\n```\n")
        f.write(scan_df.to_string(index=False))
        f.write("\n```\n\n## Patient mixing (selected K)\n\n```\n")
        f.write(mix.to_string(index=False))
        f.write("\n```\n")
        f.write("\n\n## MP top-10 signature genes\n\n")
        for mp in sorted(assign["MP"].unique(),
                         key=lambda m: int(m.replace("MP", ""))):
            top10 = (
                sigs[sigs["MP"] == mp]
                .sort_values("rank")
                .head(10)
            )
            f.write(f"### {mp} "
                    f"({(assign['MP'] == mp).sum()} GEPs, "
                    f"{int(mix.loc[mix['MP'] == mp, 'n_patients'].iloc[0])} patients)\n\n")
            f.write(", ".join(top10["gene"].tolist()) + "\n\n")
    log(f"Wrote {OUT_SUMMARY}")

    log("=" * 70)
    log("MP top-10 genes:")
    for mp in sorted(assign["MP"].unique(),
                     key=lambda m: int(m.replace("MP", ""))):
        top10 = sigs[sigs["MP"] == mp].sort_values("rank").head(10)
        n_geps = (assign["MP"] == mp).sum()
        n_pats = int(mix.loc[mix["MP"] == mp, "n_patients"].iloc[0])
        log(f"  {mp} ({n_geps} GEPs / {n_pats} pts): "
            + ", ".join(top10["gene"].tolist()))


if __name__ == "__main__":
    main()
