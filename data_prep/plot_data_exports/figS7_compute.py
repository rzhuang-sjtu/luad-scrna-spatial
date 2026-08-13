"""
Compute data for Supplementary Fig S7E (Hallmark heatmap) and S7F (GO BP).

Outputs to ${PROJECT_ROOT}/results/fig5_plot_data/:
  figS7e_hallmark_mean_by_subtype.csv      # 7 subtypes x 50 hallmarks (mean per-cell score)
  figS7e_hallmark_zscore.csv               # row-zscored variant (subtypes columns)
  figS7f_de_top200_by_subtype.csv          # subtype, gene, score, logfc, pval_adj
  figS7f_go_bp_enrichr.csv                 # subtype, term, p_adj, odds, gene_ratio, overlap, genes
"""
import os, re, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc

ANN_H5  = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad"
RAW_H5  = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad"
GMT     = "${PROJECT_ROOT}/data/gmt/MSigDB_Hallmark_2020.gmt"
OUT_DIR = "${PROJECT_ROOT}/results/fig5_plot_data"
os.makedirs(OUT_DIR, exist_ok=True)

SUBTYPE_COL = "neu_subtype"
KEEP = [
    "Neu_Inflammatory","Neu_OSM_priming","Neu_OSM_low","Neu_IFN_response",
    "Neu_Angiogenic","Neu_ECM_remodeling","Neu_Metastatic",
]

# 1) Load + normalise expression on the *raw* h5ad (full gene set)
print("[load] reading raw + annotated", flush=True)
ad_raw = sc.read_h5ad(RAW_H5)
ad_ann = sc.read_h5ad(ANN_H5)

# carry neu_subtype from annotated -> raw (same barcodes)
common = ad_raw.obs_names.intersection(ad_ann.obs_names)
print(f"[load] raw {ad_raw.shape}  annotated {ad_ann.shape}  common {len(common)}", flush=True)
ad_raw = ad_raw[common].copy()
ad_raw.obs[SUBTYPE_COL] = ad_ann.obs.loc[common, SUBTYPE_COL].values

# subset to the 7 functional subtypes
mask = ad_raw.obs[SUBTYPE_COL].isin(KEEP)
ad = ad_raw[mask].copy()
print(f"[load] kept {ad.n_obs} cells across {ad.obs[SUBTYPE_COL].nunique()} subtypes", flush=True)
print(ad.obs[SUBTYPE_COL].value_counts())

ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)

# 2) Hallmark gene-set scoring
print("[hallmark] reading gmt", flush=True)
def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 3:
                continue
            name = toks[0]
            genes = [g for g in toks[2:] if g and g != ""]
            sets[name] = genes
    return sets

hallmark = read_gmt(GMT)
print(f"[hallmark] {len(hallmark)} sets in gmt", flush=True)

scored, dropped = [], []
for name, genes in hallmark.items():
    present = [g for g in genes if g in ad.var_names]
    if len(present) < 5:
        dropped.append((name, len(present)))
        continue
    sc.tl.score_genes(ad, gene_list=present, score_name=name,
                      use_raw=False, random_state=0)
    scored.append(name)
print(f"[hallmark] scored {len(scored)} sets, dropped {len(dropped)} (<5 genes)", flush=True)

mat = ad.obs[scored + [SUBTYPE_COL]].copy()
mean_by_subtype = mat.groupby(SUBTYPE_COL, observed=True).mean().reindex(KEEP)
mean_by_subtype.index.name = "subtype"
mean_by_subtype.to_csv(os.path.join(OUT_DIR, "figS7e_hallmark_mean_by_subtype.csv"))

# row-wise z-score (per hallmark across subtypes) -> shape (subtype x hallmark)
zmat = mean_by_subtype.copy()
zmat = (zmat - zmat.mean(axis=0)) / zmat.std(axis=0).replace(0, np.nan)
zmat.to_csv(os.path.join(OUT_DIR, "figS7e_hallmark_zscore.csv"))
print("[hallmark] wrote mean + zscore CSVs", flush=True)

# 3) DE genes per subtype (vs rest)  --  for GO enrichment
print("[de] rank_genes_groups (wilcoxon, vs rest)", flush=True)
ad.obs[SUBTYPE_COL] = ad.obs[SUBTYPE_COL].astype("category")
sc.tl.rank_genes_groups(ad, groupby=SUBTYPE_COL, method="wilcoxon",
                        groups=KEEP, reference="rest", pts=True)

de_rows = []
for sub in KEEP:
    df = sc.get.rank_genes_groups_df(ad, group=sub).head(200)
    df["subtype"] = sub
    de_rows.append(df)
de_top = pd.concat(de_rows, ignore_index=True)[
    ["subtype","names","scores","logfoldchanges","pvals_adj","pct_nz_group","pct_nz_reference"]
].rename(columns={"names":"gene","logfoldchanges":"logfc","pvals_adj":"pval_adj"})
de_top.to_csv(os.path.join(OUT_DIR, "figS7f_de_top200_by_subtype.csv"), index=False)
print(f"[de] wrote {len(de_top)} rows ({de_top.subtype.nunique()} subtypes x 200)", flush=True)

# 4) GO BP enrichment via gseapy.enrichr
print("[go] enrichr GO_Biological_Process_2023", flush=True)
import gseapy

go_rows = []
for sub in KEEP:
    genes = de_top.loc[de_top.subtype == sub, "gene"].head(200).tolist()
    try:
        enr = gseapy.enrichr(
            gene_list=genes,
            gene_sets="GO_Biological_Process_2023",
            organism="human",
            outdir=None,
            no_plot=True,
        )
        df = enr.results.copy()
    except Exception as e:
        print(f"[go] {sub}: enrichr failed -- {e}", flush=True)
        continue
    df["subtype"] = sub
    df = df.rename(columns={
        "Term":"term", "P-value":"pval", "Adjusted P-value":"p_adj",
        "Old P-value":"old_p", "Old Adjusted P-value":"old_padj",
        "Odds Ratio":"odds_ratio", "Combined Score":"combined_score",
        "Genes":"genes", "Overlap":"overlap",
    })
    # gene ratio = overlap "k/N" -> k/N
    def _ratio(s):
        try:
            k, n = s.split("/")
            return int(k) / int(n)
        except Exception:
            return np.nan
    df["gene_ratio"] = df["overlap"].map(_ratio)
    df = df[["subtype","term","p_adj","pval","odds_ratio",
             "combined_score","overlap","gene_ratio","genes"]]
    go_rows.append(df.head(50))   # keep top 50 per subtype
    print(f"[go] {sub}: {len(df)} terms, top p_adj = {df['p_adj'].min():.2e}", flush=True)

go_all = pd.concat(go_rows, ignore_index=True) if go_rows else pd.DataFrame()
go_all.to_csv(os.path.join(OUT_DIR, "figS7f_go_bp_enrichr.csv"), index=False)
print(f"[go] wrote {len(go_all)} rows total", flush=True)

print("\n[done] all CSVs in", OUT_DIR)
