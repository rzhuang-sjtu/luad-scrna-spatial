"""Step 13a: GSE207422 neoadjuvant chemo-IO response (24 NSCLC pre-biopsies).

Pipeline:
  1. Load log2TPM, convert to linear TPM (2^x - 1) for ssGSEA.
  2. Score MP1-4 (top-50), EMT Hallmark, Neutrophil_core, NETs_composite.
  3. Mann-Whitney MPR vs NMPR for each score; ROC AUC.
  4. Subset analyses: all 24 NSCLC; non-squamous (Adeno+NSCLC+Carcino, n=12).

Outputs → ${WORK_ROOT}/luad_figures/fig_treatment/:
  - gse207422_mp_scores.csv         (sample × all scores + group)
  - gse207422_mpr_vs_nmpr.csv       (per score: U-test, AUC, deltas, n's)
  - gse207422_boxplot_data.csv      (long form for plotting)
  - gse207422_summary.md
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

EXPR = Path("${DATA_ROOT}/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz")
META = Path("${DATA_ROOT}/GSE207422/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
SIG = Path.home() / "luad/results/step6_mp_signatures_top100.csv"
GMT = Path.home() / "luad/data/gmt/MSigDB_Hallmark_2020.gmt"
OUT = Path("${WORK_ROOT}/luad_figures/fig_treatment")
OUT.mkdir(parents=True, exist_ok=True)

# Reuse the same gene panels as TCGA pipeline
NEUT_CORE = ["CSF3R", "FCGR3B", "CXCR1", "CXCR2", "S100A8", "S100A9",
             "MMP9", "ELANE", "CEACAM8"]
NETS_COMPOSITE = ["PADI4", "MPO", "ELANE", "CTSG", "PRTN3", "DEFA1", "DEFA3",
                  "HMGB1", "H3F3A", "DNASE1L3", "CYBB", "NCF1", "NCF2", "NCF4",
                  "S100A8", "S100A9", "S100A12", "CAMP", "LCN2", "MMP9"]
TOPN_MP = 50


def hallmark_emt() -> list[str]:
    with open(GMT) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0].lower().startswith("epithelial mesenchymal"):
                return parts[2:]
    return []


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    # ----- 1. metadata -----
    log(f"loading metadata {META}")
    meta = pd.read_excel(META)
    meta = meta[meta["Sample"].notna() & meta["Pathologic Response"].notna()].copy()
    log(f"  valid samples: {len(meta)}")
    log(f"  responses:\n{meta['Pathologic Response'].value_counts().to_string()}")

    # MPR (pCR) + MPR are responders; NMPR non-responders
    meta["response_group"] = meta["Pathologic Response"].apply(
        lambda x: "MPR" if "MPR" in str(x) and "NMPR" not in str(x) else "NMPR"
    )
    log(f"  binary groups:\n{meta['response_group'].value_counts().to_string()}")

    # ----- 2. expression -----
    log(f"loading log2TPM matrix")
    expr_log = pd.read_csv(EXPR, sep="\t", index_col=0)
    log(f"  shape: {expr_log.shape}")
    # Drop duplicates if any (collapse by max)
    if not expr_log.index.is_unique:
        expr_log = expr_log.groupby(level=0).max()
        log(f"  after dedup: {expr_log.shape}")

    # Convert log2TPM → linear TPM for ssGSEA (gseapy ssgsea uses ranks
    # so monotonic transform is fine, but keep consistent with TCGA pipeline)
    expr_lin = (np.power(2.0, expr_log) - 1.0).astype("float32")
    expr_lin = expr_lin.clip(lower=0)
    log(f"  linear TPM matrix: {expr_lin.shape}")

    # Subset to samples in metadata
    samples = [s for s in expr_lin.columns if s in set(meta["Sample"])]
    expr_lin = expr_lin[samples]
    meta = meta.set_index("Sample").loc[samples]
    log(f"  joined samples: {len(samples)}")

    # ----- 3. signatures -----
    log("preparing gene sets")
    sig_df = pd.read_csv(SIG)
    sig_df = sig_df[(sig_df["MP"].isin(["MP1","MP2","MP3","MP4"])) & (sig_df["rank"] <= TOPN_MP)]
    gene_sets = {f"MP{i}": [g for g in sig_df[sig_df["MP"]==f"MP{i}"]["gene"].tolist()
                            if g in expr_lin.index] for i in range(1,5)}
    gene_sets["EMT_Hallmark"] = [g for g in hallmark_emt() if g in expr_lin.index]
    gene_sets["Neutrophil_core"] = [g for g in NEUT_CORE if g in expr_lin.index]
    gene_sets["NETs_composite"] = [g for g in NETS_COMPOSITE if g in expr_lin.index]
    for k, v in gene_sets.items():
        log(f"  {k}: {len(v)} genes present")

    # ----- 4. ssGSEA on log2 expression (rank-based, monotonic) -----
    log("running ssGSEA")
    import gseapy as gp
    expr_use = np.log2(expr_lin + 1.0).astype("float32")  # back to log2 for stability
    ss = gp.ssgsea(
        data=expr_use, gene_sets=gene_sets, outdir=None,
        sample_norm_method="rank", no_plot=True,
        min_size=3, max_size=10000, permutation_num=0, seed=0, threads=4,
    )
    scores = ss.res2d.pivot_table(index="Name", columns="Term",
                                    values="NES").astype(float)
    scores.index.name = "Sample"
    log(f"  ssGSEA scores: {scores.shape}")

    # join with response group
    score_df = scores.join(meta[["response_group", "Pathology", "RECIST",
                                  "Residual Tumor"]], how="inner")
    score_df.to_csv(OUT / "gse207422_mp_scores.csv")

    # ----- 5. group comparison + ROC -----
    log("MPR vs NMPR comparison")
    SCORES_TO_TEST = ["MP1", "MP2", "MP3", "MP4", "EMT_Hallmark",
                       "Neutrophil_core", "NETs_composite"]
    SUBSETS = {
        "all": score_df.index.tolist(),
        "non_squamous": score_df[~score_df["Pathology"].astype(str).str.contains(
            "Squamous", na=False)].index.tolist(),
    }

    res_rows, box_rows = [], []
    for sub_name, idx in SUBSETS.items():
        sub = score_df.loc[idx]
        is_mpr = (sub["response_group"] == "MPR").astype(int).values
        n_mpr = int(is_mpr.sum())
        n_nmpr = int((1 - is_mpr).sum())
        for s in SCORES_TO_TEST:
            v = sub[s].values
            mpr_v = v[is_mpr == 1]; nmpr_v = v[is_mpr == 0]
            if len(mpr_v) < 2 or len(nmpr_v) < 2:
                continue
            U, p = mannwhitneyu(mpr_v, nmpr_v, alternative="two-sided")
            try:
                auc = roc_auc_score(is_mpr, v)
            except Exception:
                auc = np.nan
            res_rows.append({
                "subset": sub_name, "score": s,
                "n_mpr": n_mpr, "n_nmpr": n_nmpr,
                "median_mpr": float(np.median(mpr_v)),
                "median_nmpr": float(np.median(nmpr_v)),
                "delta": float(np.median(mpr_v) - np.median(nmpr_v)),
                "U": float(U), "p": float(p),
                "auc_mpr_predict": float(auc),
            })
            for samp, val, grp in zip(sub.index, v, sub["response_group"]):
                box_rows.append({"subset": sub_name, "score": s,
                                  "Sample": samp, "value": float(val),
                                  "group": grp})

    res_df = pd.DataFrame(res_rows)
    res_df.to_csv(OUT / "gse207422_mpr_vs_nmpr.csv", index=False)
    pd.DataFrame(box_rows).to_csv(OUT / "gse207422_boxplot_data.csv", index=False)

    log("\nMPR vs NMPR (all 24):")
    log(res_df[res_df["subset"]=="all"]
        .sort_values("auc_mpr_predict", ascending=False)
        .round(4).to_string(index=False))
    log("\nMPR vs NMPR (non-squamous, n={}):".format(
        len(SUBSETS["non_squamous"])))
    log(res_df[res_df["subset"]=="non_squamous"]
        .sort_values("auc_mpr_predict", ascending=False)
        .round(4).to_string(index=False))

    # ----- 6. summary md -----
    with open(OUT / "gse207422_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 13a — GSE207422 neoadjuvant chemo-IO response\n\n")
        f.write(f"- N samples (post-QC): {len(score_df)}\n")
        f.write(f"- MPR (pCR + MPR): {n_mpr}; NMPR: {n_nmpr}\n")
        f.write(f"- Pathology distribution: "
                f"{score_df['Pathology'].value_counts().to_dict()}\n\n")
        f.write("## All 24 samples (NSCLC)\n\n")
        f.write(res_df[res_df['subset']=='all']
                .sort_values('auc_mpr_predict', ascending=False)
                .round(4).to_markdown(index=False) + "\n\n")
        f.write(f"## Non-squamous subset (n={len(SUBSETS['non_squamous'])}, LUAD-enriched)\n\n")
        f.write(res_df[res_df['subset']=='non_squamous']
                .sort_values('auc_mpr_predict', ascending=False)
                .round(4).to_markdown(index=False) + "\n\n")
        f.write("## Notes\n")
        f.write("- AUC > 0.5 means score is HIGHER in MPR (responders); "
                "AUC < 0.5 means HIGHER in NMPR (non-responders).\n")
        f.write("- Sample size is small (24/12); look at effect direction "
                "and AUC rather than p-values.\n")
    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
