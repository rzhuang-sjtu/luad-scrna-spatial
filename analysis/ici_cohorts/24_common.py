"""Shared utilities for Step 24 (Fig 3H multi-cohort treatment validation)."""
from __future__ import annotations
import os, time, gzip
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.home() / "luad"
SIG_FILE = Path("${WORK_ROOT}/luad_figures/fig2/mp_signatures_top100.csv")
GMT_FILE = ROOT / "data/gmt/MSigDB_Hallmark_2020.gmt"
EXT = ROOT / "data/external"
OUT = Path("${WORK_ROOT}/luad_figures/fig_treatment")
TOPN_MP = 50

NEUT_CORE = ["CSF3R","FCGR3B","CXCR1","CXCR2","S100A8","S100A9","MMP9","ELANE","CEACAM8"]
NETS_COMPOSITE = ["PADI4","MPO","ELANE","CTSG","PRTN3","HMGB1","H3F3A","DNASE1L3",
                  "CYBB","NCF1","NCF2","NCF4","S100A8","S100A9","S100A12","CAMP","LCN2","MMP9"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def hallmark_emt() -> list[str]:
    with open(GMT_FILE) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0].lower().startswith("epithelial mesenchymal"):
                return parts[2:]
    return []


def load_mp_signatures() -> dict[str, list[str]]:
    """Top-50 genes per MP1-4 from consensus signature."""
    sig = pd.read_csv(SIG_FILE)
    sig = sig[(sig["MP"].isin(["MP1","MP2","MP3","MP4"])) & (sig["rank"] <= TOPN_MP)]
    return {f"MP{i}": sig[sig["MP"] == f"MP{i}"]["gene"].tolist() for i in range(1, 5)}


def build_gene_sets(expr_genes: pd.Index) -> tuple[dict[str, list[str]], dict[str, int]]:
    """Filter MP/aux gene sets to those present in expression matrix; return overlap count."""
    mp_sigs = load_mp_signatures()
    sets = {}
    overlaps = {}
    for k, gl in mp_sigs.items():
        present = [g for g in gl if g in expr_genes]
        sets[k] = present
        overlaps[k] = len(present)
    sets["EMT_Hallmark"] = [g for g in hallmark_emt() if g in expr_genes]
    sets["Neutrophil_core"] = [g for g in NEUT_CORE if g in expr_genes]
    sets["NETs_composite"] = [g for g in NETS_COMPOSITE if g in expr_genes]
    overlaps["EMT_Hallmark"] = len(sets["EMT_Hallmark"])
    overlaps["Neutrophil_core"] = len(sets["Neutrophil_core"])
    overlaps["NETs_composite"] = len(sets["NETs_composite"])
    return sets, overlaps


def run_ssgsea(expr_log: pd.DataFrame, gene_sets: dict[str, list[str]],
               threads: int = 4, min_size: int = 3) -> pd.DataFrame:
    """ssGSEA with rank-norm; expr_log = gene × sample log-scale matrix."""
    import gseapy as gp
    ss = gp.ssgsea(
        data=expr_log, gene_sets=gene_sets, outdir=None,
        sample_norm_method="rank", no_plot=True,
        min_size=min_size, max_size=10000, permutation_num=0,
        seed=0, threads=threads,
    )
    scores = ss.res2d.pivot_table(index="Name", columns="Term", values="NES").astype(float)
    scores.index.name = "Sample"
    return scores


def fetch_geo_soft(gse: str) -> "GEOparse.GSE":
    """GEOparse.get_GEO with cache to data/external/."""
    import GEOparse
    EXT.mkdir(parents=True, exist_ok=True)
    gse_obj = GEOparse.get_GEO(geo=gse, destdir=str(EXT), silent=True,
                                annotate_gpl=False, include_data=False)
    return gse_obj


def gsm_characteristics(gse_obj) -> pd.DataFrame:
    """Build a DataFrame: GSM × {title, characteristics_ch1 list joined}."""
    rows = []
    for gsm_name, gsm in gse_obj.gsms.items():
        d = {"GSM": gsm_name,
             "title": (gsm.metadata.get("title") or [""])[0]}
        chars = gsm.metadata.get("characteristics_ch1", [])
        for c in chars:
            if ":" in c:
                k, v = c.split(":", 1)
                k = k.strip().lower().replace(" ", "_")
                v = v.strip()
                d[k] = v
        rows.append(d)
    return pd.DataFrame(rows)


def wilcoxon_auc(values: np.ndarray, labels_pos: np.ndarray) -> dict:
    """Mann-Whitney + ROC AUC. labels_pos: 1 for "responder" (positive class)."""
    from scipy.stats import mannwhitneyu
    from sklearn.metrics import roc_auc_score
    pos = values[labels_pos == 1]
    neg = values[labels_pos == 0]
    if len(pos) < 2 or len(neg) < 2:
        return dict(n_pos=len(pos), n_neg=len(neg), median_pos=np.nan,
                    median_neg=np.nan, delta=np.nan, U=np.nan, p=np.nan, auc=np.nan)
    U, p = mannwhitneyu(pos, neg, alternative="two-sided")
    try:
        auc = roc_auc_score(labels_pos, values)
    except Exception:
        auc = np.nan
    return dict(n_pos=int(len(pos)), n_neg=int(len(neg)),
                median_pos=float(np.median(pos)), median_neg=float(np.median(neg)),
                delta=float(np.median(pos) - np.median(neg)),
                U=float(U), p=float(p), auc=float(auc))


def write_outputs(gse: str, score_df: pd.DataFrame, group_col: str,
                  pos_label: str, neg_label: str,
                  scores_to_test: list[str] | None = None,
                  extra_meta_cols: list[str] | None = None):
    """Write standard 3 files: mp_scores, response_comparison, boxplot_data."""
    OUT.mkdir(parents=True, exist_ok=True)
    if scores_to_test is None:
        scores_to_test = ["MP1", "MP2", "MP3", "MP4", "EMT_Hallmark",
                          "Neutrophil_core", "NETs_composite"]

    score_df.to_csv(OUT / f"{gse}_mp_scores.csv")

    # Build comparison
    is_pos = (score_df[group_col].astype(str) == pos_label).astype(int).values
    rows = []
    for s in scores_to_test:
        if s not in score_df.columns:
            continue
        v = score_df[s].values.astype(float)
        stats = wilcoxon_auc(v, is_pos)
        stats["score"] = s
        stats["pos_group"] = pos_label
        stats["neg_group"] = neg_label
        rows.append(stats)
    comp = pd.DataFrame(rows)[["score", "pos_group", "neg_group", "n_pos", "n_neg",
                                 "median_pos", "median_neg", "delta",
                                 "U", "p", "auc"]]
    comp.to_csv(OUT / f"{gse}_response_comparison.csv", index=False)

    # Long-format boxplot data
    box_rows = []
    for s in scores_to_test:
        if s not in score_df.columns:
            continue
        for samp, val, grp in zip(score_df.index,
                                    score_df[s].values, score_df[group_col]):
            box_rows.append({"Sample": samp, "score": s,
                              "value": float(val), "group": str(grp)})
    pd.DataFrame(box_rows).to_csv(OUT / f"{gse}_boxplot_data.csv", index=False)
    return comp
