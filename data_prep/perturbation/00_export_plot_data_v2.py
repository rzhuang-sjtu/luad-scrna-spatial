"""Export ALL plot-ready CSVs needed by Fig 8 v2 (500-cell, leads = SEC61G/SRSF9/ANGPTL4).

Output dir: ${PROJECT_ROOT}/results/fig8_plot_data/v2_500/
The R panel scripts read from this dir only.
"""
from pathlib import Path
import re
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
warnings.filterwarnings("ignore")

OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/v2_500")
OUT.mkdir(parents=True, exist_ok=True)

PERT_500 = Path("${PROJECT_ROOT}/results/fig8_geneformer/perturb_500")
DEPMAP   = Path("${DATA_ROOT}/depmap/24Q2")
TCGA     = Path("${DATA_ROOT}/TCGA_LUAD_analysis")
ROI_E    = Path("${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad")
ROI_O    = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/cohort_with_roi.h5ad")
PD_      = Path("${WORK_ROOT}/luad_figures/fig_treatment")

LEADS = ["SEC61G", "SRSF9", "ANGPTL4"]
ENSG  = {"ANGPTL4": "ENSG00000167772", "EGLN3": "ENSG00000129521",
         "SEC61G": "ENSG00000132432", "SLC2A1": "ENSG00000117394",
         "SRSF9":  "ENSG00000111786"}
TRANSITIONS = ["macro_spp1_to_c1qc", "mal_mp3_to_mp1", "neu_osm_priming_to_low"]
N_MIN, FDR = 5, 0.05

print("[8A] venn sets")
def load_filtered(t):
    df = pd.read_csv(PERT_500 / t / f"{t}_stats.csv", index_col=0)
    df = df[(df["Sig"] == 1) & (df["Shift_to_goal_end"] > 0)
            & (df["Goal_end_FDR"] < FDR) & (df["N_Detections"] >= N_MIN)].copy()
    return df.sort_values("Shift_to_goal_end", ascending=False).reset_index(drop=True)
filt = {t: load_filtered(t) for t in TRANSITIONS}
sizes = pd.DataFrame([{"transition": t, "n_filtered": len(filt[t])} for t in TRANSITIONS])
sizes.to_csv(OUT / "8A_pool_sizes.csv", index=False)

venn_rows = []
for n in (50, 100, 200, 300):
    sets = {t: set(filt[t].head(n)["Gene_name"]) for t in TRANSITIONS}
    a, b, c = sets["macro_spp1_to_c1qc"], sets["mal_mp3_to_mp1"], sets["neu_osm_priming_to_low"]
    venn_rows.append({"top_N": n,
                      "only_macro": len(a - b - c), "only_mal": len(b - a - c), "only_neu": len(c - a - b),
                      "macro_AND_mal_only":  len((a & b) - c),
                      "macro_AND_neu_only":  len((a & c) - b),
                      "mal_AND_neu_only":    len((b & c) - a),
                      "all_three":           len(a & b & c),
                      "abc_genes": ";".join(sorted(a & b & c)),
                      "size_macro": len(a), "size_mal": len(b), "size_neu": len(c)})
pd.DataFrame(venn_rows).to_csv(OUT / "8A_venn_subsets.csv", index=False)

print("[8B-C] DepMap")
model = pd.read_csv(DEPMAP / "Model.csv", low_memory=False)
luad_models = set(model.loc[model["OncotreeCode"] == "LUAD", "ModelID"])

crispr_header = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", nrows=0).columns.tolist()
def name_of(c):
    m = re.match(r"^([A-Za-z0-9\-]+)\s*\(\d+\)$", c)
    return m.group(1) if m else c
sym2col = {name_of(c): c for c in crispr_header[1:]}
sel = [crispr_header[0]] + [sym2col[g] for g in LEADS if g in sym2col]
crispr = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", usecols=sel)
crispr = crispr.rename(columns={crispr_header[0]: "ModelID"})
crispr = crispr.rename(columns={sym2col[g]: g for g in LEADS if g in sym2col})
crispr["is_LUAD"] = crispr["ModelID"].isin(luad_models)
crispr["group"]   = np.where(crispr["is_LUAD"], "LUAD", "non-LUAD")

dep_long = crispr.melt(id_vars=["ModelID", "is_LUAD", "group"],
                       value_vars=LEADS, var_name="gene", value_name="gene_effect").dropna()
dep_long.to_csv(OUT / "8B_depmap_long.csv", index=False)

dep_stats = []
for g in LEADS:
    a = crispr.loc[crispr["is_LUAD"], g].dropna()
    b = crispr.loc[~crispr["is_LUAD"], g].dropna()
    p = stats.mannwhitneyu(a, b, alternative="less").pvalue if len(a) >= 5 and len(b) >= 5 else np.nan
    dep_stats.append({"gene": g, "luad_mean": a.mean(), "luad_n": len(a),
                      "other_mean": b.mean(), "other_n": len(b),
                      "delta_LUAD_minus_other": a.mean() - b.mean(),
                      "mw_p_LUAD_lt_other": p})
pd.DataFrame(dep_stats).to_csv(OUT / "8B_depmap_stats.csv", index=False)

# DepMap expression (TPM logp1) for 8G scatter — file has meta cols (SequencingID,
# ModelConditionID, ModelID, IsDefaultEntryForMC, IsDefaultEntryForModel) before genes
print("  loading DepMap expression for 8G scatter ...")
expr_header = pd.read_csv(DEPMAP / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                          nrows=0).columns.tolist()
META = {"SequencingID", "ModelConditionID", "ModelID",
        "IsDefaultEntryForMC", "IsDefaultEntryForModel"}
expr_id = expr_header[0]  # unnamed first col
e_sym2col = {name_of(c): c for c in expr_header
             if c not in META and c != expr_id}
e_use = [expr_id, "ModelID", "IsDefaultEntryForModel"] + \
        [e_sym2col[g] for g in LEADS if g in e_sym2col]
expr = pd.read_csv(DEPMAP / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                   usecols=e_use)
expr = expr[expr["IsDefaultEntryForModel"] == "Yes"].copy()
expr = expr.rename(columns={e_sym2col[g]: g for g in LEADS if g in e_sym2col})
expr_long = expr.melt(id_vars=["ModelID"], value_vars=LEADS,
                      var_name="gene", value_name="log2_TPM_p1").dropna()
expr_long["ModelID"] = expr_long["ModelID"].astype(str)
expr_long["is_LUAD"] = expr_long["ModelID"].isin(luad_models)

# join with effect
dep_long2 = dep_long.copy(); dep_long2["ModelID"] = dep_long2["ModelID"].astype(str)
g8 = pd.merge(expr_long, dep_long2[["ModelID", "gene", "gene_effect"]],
              on=["ModelID", "gene"], how="inner")
print(f"  8G rows after merge: {len(g8)} (LUAD={g8['is_LUAD'].sum()})")
g8.to_csv(OUT / "8G_expr_vs_effect.csv", index=False)

print("[8D-F] TCGA T vs N")
cln = pd.read_csv(TCGA / "TCGA_LUAD_clinical.csv").rename(columns={"sample_barcode": "sample"})
tpm_rows = []
for chunk in pd.read_csv(TCGA / "TCGA_LUAD_TPM_matrix.csv", chunksize=5000):
    chunk = chunk.rename(columns={chunk.columns[0]: "gene"})
    sub = chunk[chunk["gene"].isin(LEADS)]
    if len(sub): tpm_rows.append(sub)
tpm = pd.concat(tpm_rows, ignore_index=True).set_index("gene")

st_map = dict(zip(cln["sample"], cln["sample_type"]))
rows = []
for g in tpm.index:
    for s, v in tpm.loc[g].items():
        st = st_map.get(s, "Unknown")
        if st in ("Primary Tumor", "Solid Tissue Normal"):
            rows.append({"gene": g, "sample": s, "log2_TPM_p1": np.log2(v + 1),
                         "type": "Tumor" if st == "Primary Tumor" else "Normal"})
tn_long = pd.DataFrame(rows)
tn_long.to_csv(OUT / "8D_tcga_TvN_long.csv", index=False)

tn_stats = []
for g in LEADS:
    sub = tn_long[tn_long["gene"] == g]
    t = sub.loc[sub["type"] == "Tumor", "log2_TPM_p1"]
    n = sub.loc[sub["type"] == "Normal", "log2_TPM_p1"]
    u, p = stats.mannwhitneyu(t, n, alternative="two-sided")
    tn_stats.append({"gene": g, "tumor_mean": t.mean(), "tumor_n": len(t),
                     "normal_mean": n.mean(), "normal_n": len(n),
                     "log2FC_T_minus_N": t.mean() - n.mean(), "wilcoxon_p": p})
pd.DataFrame(tn_stats).to_csv(OUT / "8D_tcga_TvN_stats.csv", index=False)

print("[8O/P] KM survival")
tumor_samples = [s for s in cln.loc[cln["sample_type"] == "Primary Tumor", "sample"]
                 if s in tpm.columns]
surv = cln[cln["sample"].isin(tumor_samples)].copy()
surv["event"] = (surv["vital_status"] == "Dead").astype(int)
surv["time"]  = np.where(surv["event"] == 1, surv["days_to_death"], surv["days_to_last_follow_up"])
surv = surv.dropna(subset=["time"]); surv = surv[surv["time"] > 0].copy()

km_rows = []
for g in ("SRSF9", "SEC61G"):
    if g not in tpm.index: continue
    expr_g = tpm.loc[g, surv["sample"].values].values
    s = surv.copy()
    s["expr"] = np.log2(expr_g + 1)
    med = s["expr"].median()
    s["group"] = np.where(s["expr"] >= med, "High", "Low")
    s["gene"] = g
    km_rows.append(s[["gene", "sample", "time", "event", "expr", "group"]])
    high = s[s["group"] == "High"]; low = s[s["group"] == "Low"]
    lr = logrank_test(high["time"], low["time"], high["event"], low["event"])
    cph = CoxPHFitter()
    cph.fit(s[["time","event","expr"]], duration_col="time", event_col="event")
    sm = cph.summary.loc["expr"]
    pd.DataFrame([{"gene": g, "n_high": len(high), "n_low": len(low),
                   "logrank_p": lr.p_value,
                   "HR": np.exp(sm["coef"]),
                   "HR_low": np.exp(sm["coef lower 95%"]),
                   "HR_high": np.exp(sm["coef upper 95%"]),
                   "cox_p": sm["p"]}]).to_csv(OUT / f"8{'O' if g=='SRSF9' else 'P'}_km_{g}_stats.csv", index=False)
pd.concat(km_rows).to_csv(OUT / "8OP_km_long.csv", index=False)

print("[8H/I-L/S11ABCIJ + atlases] ST E-MTAB-13530 + Okamura spot-level")
ad = sc.read_h5ad(ROI_E)
print(f"  E-MTAB loaded: {ad.shape}")

# section-level mean MP scores → R/NR surrogate (E-MTAB)
sec_mp = ad.obs.groupby("sample")[["MP3_score", "MP4_score"]].mean().reset_index()
sec_mp.to_csv(OUT / "8I_section_MP_scores.csv", index=False)
sec_R  = sec_mp.sort_values("MP4_score", ascending=False).head(2)["sample"].tolist()
sec_NR = sec_mp.sort_values("MP3_score", ascending=False).head(2)["sample"].tolist()
used = set(sec_R) | set(sec_NR)
sec_extra = [s for s in sec_mp["sample"] if s not in used][:3]
print(f"  E-MTAB R-surrogate (8I/J): {sec_R}")
print(f"  E-MTAB NR-surrogate (8K/L): {sec_NR}")
print(f"  E-MTAB S11A-C extra: {sec_extra}")
all_emtab_sections = sorted(sec_mp["sample"].astype(str).unique().tolist())

# Okamura cohort
ad_o = sc.read_h5ad(ROI_O)
print(f"  Okamura loaded: {ad_o.shape}")
sec_mp_o = ad_o.obs.groupby("sample")[["MP3_score", "MP4_score"]].mean().reset_index()
sec_mp_o.to_csv(OUT / "okamura_section_MP_scores.csv", index=False)
sec_R_o  = sec_mp_o.sort_values("MP4_score", ascending=False).head(1)["sample"].tolist()
sec_NR_o = sec_mp_o.sort_values("MP3_score", ascending=False).head(1)["sample"].tolist()
print(f"  Okamura R-surrogate (S11I): {sec_R_o}")
print(f"  Okamura NR-surrogate (S11J): {sec_NR_o}")
all_okamura_sections = sorted(sec_mp_o["sample"].astype(str).unique().tolist())

# panel assignments for INDIVIDUAL representative panels (1 per section, 1x3 gene strip)
panel_assignments = (
    [(s, "E-MTAB-13530", "8I", "R-surrogate")  for s in sec_R[:1]] +
    [(s, "E-MTAB-13530", "8J", "R-surrogate")  for s in sec_R[1:2]] +
    [(s, "E-MTAB-13530", "8K", "NR-surrogate") for s in sec_NR[:1]] +
    [(s, "E-MTAB-13530", "8L", "NR-surrogate") for s in sec_NR[1:2]] +
    [(s, "E-MTAB-13530", f"S11{chr(65+i)}", "extra") for i, s in enumerate(sec_extra)] +
    [(s, "Okamura",      "S11I", "R-surrogate")  for s in sec_R_o[:1]] +
    [(s, "Okamura",      "S11J", "NR-surrogate") for s in sec_NR_o[:1]]
)
pd.DataFrame(panel_assignments, columns=["sample","cohort","panel","kind"]
            ).to_csv(OUT / "8I_panel_assignments.csv", index=False)

# spot-level data for ALL sections (28) × 3 leads. R picks the subset it needs.
# Output raw full-res spatial coords; R applies per-section tissue_hires_scalef.
def _dump_cohort(adobj, cohort_name, sections):
    rows = []
    for s in sections:
        mask = adobj.obs["sample"].astype(str) == s
        sub = adobj[mask]
        if sub.n_obs == 0:
            continue
        if "spatial" in sub.obsm:
            xy = np.asarray(sub.obsm["spatial"])
        else:
            xy = sub.obs[["x", "y"]].values
        for g in LEADS:
            if g not in sub.var_names:
                continue
            e = sub[:, g].X
            e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
            for k in range(len(e)):
                rows.append({"cohort": cohort_name, "sample": s, "gene": g,
                             "spatial1": float(xy[k, 0]),
                             "spatial2": float(xy[k, 1]),
                             "expr": float(e[k])})
    return rows

spot_rows = (_dump_cohort(ad,   "E-MTAB-13530", all_emtab_sections) +
             _dump_cohort(ad_o, "Okamura",      all_okamura_sections))
pd.DataFrame(spot_rows).to_csv(OUT / "8I_spot_long.csv", index=False)
print(f"  spot CSV rows: {len(spot_rows)}  (E-MTAB sections: {len(all_emtab_sections)}, "
      f"Okamura sections: {len(all_okamura_sections)})")

# 8H: co-expression on E-MTAB whole cohort (UMAP-ish layout from 2D PCA on log-expr or
#     just X_umap if present). Inspect.
print("  8H co-expression representation ...")
emb = None
for k in ("X_umap", "X_umap_harmony", "X_pca"):
    if k in ad.obsm: emb = ad.obsm[k]; emb_name = k; break
if emb is None:
    # fall back: 2D PCA on lead-gene log-expr
    Xg = np.column_stack([
        (ad[:, g].X.toarray().flatten() if hasattr(ad[:, g].X, "toarray") else
         np.asarray(ad[:, g].X).flatten())
        for g in LEADS if g in ad.var_names
    ])
    if Xg.shape[1] >= 2:
        from numpy.linalg import svd
        m = Xg - Xg.mean(0)
        u, s_, vt = svd(m, full_matrices=False)
        emb = u[:, :2] * s_[:2]
        emb_name = "leadPCA"
print(f"  embedding: {emb_name}")
emb2 = emb[:, :2]
co_rows = []
for g in LEADS:
    if g not in ad.var_names: continue
    e = ad[:, g].X
    e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
    n = len(e)
    co_rows.append(pd.DataFrame({"gene": [g] * n, "x": emb2[:, 0], "y": emb2[:, 1],
                                 "expr": e, "sample": ad.obs["sample"].values}))
pd.concat(co_rows).to_csv(OUT / "8H_coexpr_embedding.csv", index=False)
pd.DataFrame([{"embedding": emb_name}]).to_csv(OUT / "8H_embedding_info.csv", index=False)

print("[8M/N/S11D-F] treatment cohorts")
def load_expr(path, kind):
    if kind == "log2TPM_symbol":
        return pd.read_csv(path, sep="\t", index_col=0)
    elif kind == "counts_symbol":
        df = pd.read_csv(path, sep="\t", index_col=0)
        if not df.index.is_unique: df = df.groupby(level=0).max()
        cpm = df.div(df.sum(axis=0), axis=1) * 1e6
        return np.log2(cpm + 1)
    elif kind == "TPM_ENSG":
        df = pd.read_csv(path, sep="\t", index_col=0)
        df.index = df.index.astype(str).str.split(".").str[0]
        if not df.index.is_unique: df = df.groupby(level=0).max()
        return np.log2(df + 1)

# C1 = 8M (GSE207422 MPR vs NMPR — 3 lead boxes)
C1 = dict(name="GSE207422", expr="${DATA_ROOT}/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
          kind="log2TPM_symbol", scores=PD_/"gse207422_mp_scores.csv", pos="MPR", neg="NMPR")
e1 = load_expr(C1["expr"], C1["kind"])
sc1 = pd.read_csv(C1["scores"]).set_index("Sample")
samples = [s for s in e1.columns if s in sc1.index]
e1 = e1[samples]; resp = sc1.loc[samples, "response_group"]
keep = resp.isin([C1["pos"], C1["neg"]]); e1 = e1.loc[:, keep.values]; resp = resp[keep]
rows = []
for g in LEADS:
    if g in e1.index:
        for s in e1.columns:
            rows.append({"gene": g, "sample": s, "expr": float(e1.loc[g, s]),
                         "response": resp[s]})
pd.DataFrame(rows).to_csv(OUT / "8M_GSE207422_long.csv", index=False)

stat_rows = []
for g in LEADS:
    sub = pd.DataFrame(rows); sub = sub[sub["gene"] == g]
    a = sub.loc[sub["response"] == C1["pos"], "expr"]
    b = sub.loc[sub["response"] == C1["neg"], "expr"]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan)
    auc = (u / (len(a) * len(b))) if not np.isnan(p) else np.nan
    stat_rows.append({"gene": g, "n_pos": len(a), "n_neg": len(b),
                      "log2FC_pos_minus_neg": a.mean() - b.mean(),
                      "p": p, "auc_pos_vs_neg": auc})
pd.DataFrame(stat_rows).to_csv(OUT / "8M_GSE207422_stats.csv", index=False)

# C2 = 8N (GSE126044 R vs NR full-genome volcano + 3 leads highlighted)
C2 = dict(name="GSE126044", expr="${DATA_ROOT}/GSE126044/GSE126044_counts.txt.gz",
          kind="counts_symbol", scores=PD_/"gse126044_mp_scores.csv", pos="R", neg="NR")
e2 = load_expr(C2["expr"], C2["kind"])
sc2 = pd.read_csv(C2["scores"]).set_index("Sample")
samples = [s for s in e2.columns if s in sc2.index]
e2 = e2[samples]; resp = sc2.loc[samples, "response_group"]
keep = resp.isin([C2["pos"], C2["neg"]]); e2 = e2.loc[:, keep.values]; resp = resp[keep]
mask_pos = (resp.values == C2["pos"]); mask_neg = (resp.values == C2["neg"])

vol = []
for g in e2.index:
    a = e2.loc[g, mask_pos].astype(float).values
    b = e2.loc[g, mask_neg].astype(float).values
    if (a == a[0]).all() and (b == b[0]).all() and a[0] == b[0]: continue
    try:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    except Exception:
        p = np.nan
    vol.append({"gene": g, "log2FC_R_minus_NR": float(a.mean() - b.mean()), "p": float(p) if not np.isnan(p) else np.nan})
vol_df = pd.DataFrame(vol).dropna(subset=["p"])
vol_df["nlog10p"] = -np.log10(vol_df["p"])
vol_df["is_lead"] = vol_df["gene"].isin(LEADS)
vol_df.to_csv(OUT / "8N_GSE126044_volcano.csv", index=False)

# C3 = S11D-F (GSE135222 R vs NR — 3 lead boxes)
C3 = dict(name="GSE135222", expr="${DATA_ROOT}/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
          kind="TPM_ENSG", scores=PD_/"gse135222_mp_scores.csv", pos="R", neg="NR")
e3 = load_expr(C3["expr"], C3["kind"])
sc3 = pd.read_csv(C3["scores"]).set_index("Sample")
samples = [s for s in e3.columns if s in sc3.index]
e3 = e3[samples]; resp = sc3.loc[samples, "response_group"]
keep = resp.isin([C3["pos"], C3["neg"]]); e3 = e3.loc[:, keep.values]; resp = resp[keep]
rows3, stat3 = [], []
for g in LEADS:
    eg = ENSG.get(g)
    if eg in e3.index:
        for s in e3.columns:
            rows3.append({"gene": g, "sample": s, "expr": float(e3.loc[eg, s]),
                          "response": resp[s]})
pd.DataFrame(rows3).to_csv(OUT / "S11D_GSE135222_long.csv", index=False)
for g in LEADS:
    sub = pd.DataFrame(rows3); sub = sub[sub["gene"] == g]
    if not len(sub):
        stat3.append({"gene": g, "n_pos": 0, "n_neg": 0,
                      "log2FC_pos_minus_neg": np.nan, "p": np.nan, "auc_pos_vs_neg": np.nan}); continue
    a = sub.loc[sub["response"] == C3["pos"], "expr"]
    b = sub.loc[sub["response"] == C3["neg"], "expr"]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan)
    auc = (u / (len(a) * len(b))) if not np.isnan(p) else np.nan
    stat3.append({"gene": g, "n_pos": len(a), "n_neg": len(b),
                  "log2FC_pos_minus_neg": a.mean() - b.mean(),
                  "p": p, "auc_pos_vs_neg": auc})
pd.DataFrame(stat3).to_csv(OUT / "S11D_GSE135222_stats.csv", index=False)

print(f"\nALL CSV exports done. Files in: {OUT}")
print("\n".join(f"  {p.name}" for p in sorted(OUT.iterdir()) if p.suffix == ".csv"))
