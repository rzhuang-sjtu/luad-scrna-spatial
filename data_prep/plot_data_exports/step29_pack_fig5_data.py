"""step29: package all Figure 5 panel plotting data into R-ready CSVs.

Output dir: ~/luad/results/fig5_plot_data/
Format: long format where appropriate (ggplot2 friendly), wide where natural (heatmap matrix).
All CSVs use functional Neu_* labels.
"""
import os, time, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import spearmanr
from scipy.sparse import issparse

t0 = time.time()
DATA = Path("${PROJECT_ROOT}/data/processed")
RES = Path("${PROJECT_ROOT}/results")
OUT = RES / "fig5_plot_data"
OUT.mkdir(exist_ok=True, parents=True)

# Functional rename map (Salcher → functional)
RENAME = {
    "TAN-1": "Neu_Inflammatory",
    "TAN-2": "Neu_IFN_response",
    "TAN-3": "Neu_Angiogenic",
    "TAN-4": "Neu_Metastatic",
    "NAN-1": "Neu_ECM_remodeling",
    "NAN-2": "Neu_OSM_priming",
    "NAN-3": "Neu_OSM_low",
    "Neutrophils": "Neu_unclassified",
}
ORDER = ["Neu_Inflammatory", "Neu_Angiogenic", "Neu_Metastatic",
         "Neu_IFN_response", "Neu_OSM_priming", "Neu_OSM_low",
         "Neu_ECM_remodeling", "Neu_unclassified"]

def assert_clean(df, label):
    nas = [c for c in df.columns if pd.isna(c) or str(c).strip() == ""]
    if nas:
        raise ValueError(f"{label}: NA/empty column names: {nas}")
    if df.columns.duplicated().any():
        dup = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(f"{label}: duplicate column names: {dup}")

# 5A: UMAP metadata
print("[5A] umap metadata")
ann = sc.read_h5ad(DATA / "luad_neutrophil_own_annotated.h5ad")
um = ann.obsm["X_umap"]
df_a = pd.DataFrame({
    "barcode":      ann.obs.index.astype(str),
    "UMAP1":        um[:, 0],
    "UMAP2":        um[:, 1],
    "neu_subtype":  ann.obs["neu_subtype"].astype(str).values,
    "scanvi_predicted": ann.obs["scanvi_predicted"].astype(str).values,  # archeology
    "scanvi_uncertainty": ann.obs["scanvi_uncertainty"].astype(float).values,
    "leiden_0.6":   ann.obs["leiden_0.6"].astype(str).values,
    "tissue_type":  ann.obs["tissue_type"].astype(str).values,
    "dataset":      ann.obs["dataset"].astype(str).values,
    "patient_id":   ann.obs["patient_id"].astype(str).values,
    "sample_id":    ann.obs["sample_id"].astype(str).values,
})
# rename leiden_0.6 column key (dot in colname is annoying for R)
df_a = df_a.rename(columns={"leiden_0.6": "leiden_0_6"})
assert_clean(df_a, "5A")
df_a.to_csv(OUT / "fig5a_umap_metadata.csv.gz", index=False, compression="gzip")
print(f"  {df_a.shape} → fig5a_umap_metadata.csv.gz")

# 5B: tissue proportion (long format, normalized to 100% per tissue)
print("\n[5B] tissue proportion")
xt = pd.crosstab(ann.obs["neu_subtype"].astype(str),
                 ann.obs["tissue_type"].astype(str),
                 normalize="columns") * 100
xt = xt.reindex(ORDER).fillna(0)
df_b = xt.reset_index().melt(id_vars="neu_subtype", var_name="tissue_type",
                              value_name="proportion_pct")
# also include raw counts for ggplot label use
xt_n = pd.crosstab(ann.obs["neu_subtype"].astype(str),
                   ann.obs["tissue_type"].astype(str)).reindex(ORDER).fillna(0).astype(int)
df_b_n = xt_n.reset_index().melt(id_vars="neu_subtype", var_name="tissue_type",
                                  value_name="n_cells")
df_b = df_b.merge(df_b_n, on=["neu_subtype", "tissue_type"], how="left")
# add total per tissue
tissue_totals = ann.obs["tissue_type"].value_counts().rename("tissue_total_cells").to_frame().reset_index()
tissue_totals.columns = ["tissue_type", "tissue_total_cells"]
df_b = df_b.merge(tissue_totals, on="tissue_type", how="left")
assert_clean(df_b, "5B")
df_b.to_csv(OUT / "fig5b_tissue_proportion.csv", index=False)
print(f"  {df_b.shape} → fig5b_tissue_proportion.csv")

# helper: long-format dot plot (mean_expression + pct_expressing per subtype × gene)
def make_dot_table(adata, gene_list, group_col):
    """Compute mean log-normalized expression and % expressing per (group, gene)."""
    genes = [g for g in gene_list if g in adata.var_names]
    if not genes:
        return pd.DataFrame()
    X = adata[:, genes].X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(X, columns=genes, index=adata.obs.index)
    df[group_col] = adata.obs[group_col].astype(str).values
    means = df.groupby(group_col, observed=True)[genes].mean()
    fracs = df.groupby(group_col, observed=True)[genes].apply(lambda d: (d > 0).mean())
    out = []
    for g in genes:
        for grp in means.index:
            out.append({
                group_col:     grp,
                "gene":        g,
                "mean_expression":  float(means.loc[grp, g]),
                "pct_expressing":   float(fracs.loc[grp, g]) * 100.0,
            })
    return pd.DataFrame(out)

# 5C: canonical markers per subtype (top 5-8 each, function-aligned)
print("\n[5C] canonical markers dot")
# load raw 9698 for accurate expression; merge labels
raw = sc.read_h5ad(DATA / "luad_neutrophil_own_raw.h5ad")
raw.obs["neu_subtype"] = ann.obs.loc[raw.obs.index, "neu_subtype"].astype(str).values
sc.pp.normalize_total(raw, target_sum=1e4)
sc.pp.log1p(raw)

# Per-subtype canonical markers (literature, filtered to what's in 9698 gene set)
CANON = {
    "Neu_Inflammatory":   ["CXCL8", "IL1B", "IL1A", "PLAU", "PLAUR", "BCL2A1", "ICAM1", "VEGFA"],
    "Neu_Angiogenic":     ["VEGFA", "CXCL2", "CXCL8", "PTGS2", "G0S2", "CXCR4"],
    "Neu_Metastatic":     ["S100A8", "S100A9", "LCN2", "CCL3", "CXCL2", "AREG"],
    "Neu_IFN_response":   ["IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "RSAD2", "OAS1", "IFI6"],
    "Neu_OSM_priming":    ["OSM", "PLAUR", "VCAN", "CD55", "MNDA", "CXCL8"],
    "Neu_OSM_low":        ["OSM", "LCN2", "MMP9", "RETN"],
    "Neu_ECM_remodeling": ["MMP9", "LCN2", "CSF3R", "FPR1", "FCGR2A"],
    "Neu_unclassified":   ["AREG", "EREG", "FN1", "CXCR4"],
}
canon_genes_all = sorted({g for v in CANON.values() for g in v})
df_c = make_dot_table(raw, canon_genes_all, "neu_subtype")
# add 'category' (which subtype this gene is canonical for) — first hit only
gene_to_cat = {}
for cat, genes in CANON.items():
    for g in genes:
        if g not in gene_to_cat:
            gene_to_cat[g] = cat
df_c["canonical_for"] = df_c["gene"].map(gene_to_cat)
# enforce subtype display order
df_c["neu_subtype"] = pd.Categorical(df_c["neu_subtype"], categories=ORDER, ordered=True)
df_c = df_c.sort_values(["neu_subtype", "gene"]).reset_index(drop=True)
assert_clean(df_c, "5C")
df_c.to_csv(OUT / "fig5c_canonical_markers.csv", index=False)
print(f"  {df_c.shape} → fig5c_canonical_markers.csv ; n_genes={df_c['gene'].nunique()}")

# 5D: EMT ligand dot (25 ligands)
print("\n[5D] EMT ligand dot")
EMT_LIGANDS_ORDERED = [
    # TGFb
    "TGFB1",
    # TNF/IL axis
    "TNF", "IL6", "IL1B", "IL1A", "OSM",
    # chemokines
    "CXCL8", "CXCL1", "CXCL2", "CCL2", "CCL3", "CCL4", "CCL5",
    # MMP
    "MMP9",
    # angio/EMT
    "VEGFA", "VEGFB", "FN1", "SPP1", "SERPINE1", "PLAU", "PLAUR",
    # other
    "AREG", "EREG", "WNT5A", "PDGFB",
]
df_d = make_dot_table(raw, EMT_LIGANDS_ORDERED, "neu_subtype")
# add gene_family for facet/grouping
GENE_FAMILY = {
    "TGFB1": "TGFb",
    "TNF": "TNF/IL", "IL6": "TNF/IL", "IL1B": "TNF/IL", "IL1A": "TNF/IL", "OSM": "TNF/IL",
    "CXCL8": "Chemokine", "CXCL1": "Chemokine", "CXCL2": "Chemokine",
    "CCL2": "Chemokine", "CCL3": "Chemokine", "CCL4": "Chemokine", "CCL5": "Chemokine",
    "MMP9": "MMP",
    "VEGFA": "Angio_EMT", "VEGFB": "Angio_EMT", "FN1": "Angio_EMT", "SPP1": "Angio_EMT",
    "SERPINE1": "Angio_EMT", "PLAU": "Angio_EMT", "PLAUR": "Angio_EMT",
    "AREG": "EGF_other", "EREG": "EGF_other", "WNT5A": "EGF_other", "PDGFB": "EGF_other",
}
df_d["gene_family"] = df_d["gene"].map(GENE_FAMILY)
df_d["neu_subtype"] = pd.Categorical(df_d["neu_subtype"], categories=ORDER, ordered=True)
# preserve gene panel order
df_d["gene"] = pd.Categorical(df_d["gene"], categories=EMT_LIGANDS_ORDERED, ordered=True)
df_d = df_d.sort_values(["neu_subtype", "gene"]).reset_index(drop=True)
assert_clean(df_d, "5D")
df_d.to_csv(OUT / "fig5d_emt_ligand_dotplot.csv", index=False)
print(f"  {df_d.shape} → fig5d_emt_ligand_dotplot.csv ; n_genes={df_d['gene'].nunique()}")

# 5E: EMT logFC matrix (gene × subtype, value = logFC vs rest)
print("\n[5E] EMT logFC matrix")
src = pd.read_csv(RES / "step25f_emt_logfc_matrix.csv", index_col=0)
# rename columns from TAN/NAN → functional
src.columns = [RENAME.get(c, c) for c in src.columns]
# enforce column order, fill missing with 0
ord_cols_present = [c for c in ORDER if c in src.columns]
src = src[ord_cols_present]
# drop "Neu_unclassified" if user wants only 7 functional types (keep though for completeness)
src = src.reset_index().rename(columns={"index": "gene"})
assert_clean(src, "5E")
src.to_csv(OUT / "fig5e_emt_logfc_matrix.csv", index=False)
print(f"  {src.shape} → fig5e_emt_logfc_matrix.csv (wide: gene × subtype)")

# also long form for convenience
src_long = src.melt(id_vars="gene", var_name="neu_subtype", value_name="logFC")
src_long.to_csv(OUT / "fig5e_emt_logfc_long.csv", index=False)

# 5F: per-patient single-cell Neu × MP correlation (rho + p + padj)
print("\n[5F] single-cell Neu × MP correlation")
sc_corr = pd.read_csv(RES / "step26_tan_mp_spearman.csv")
# original col tan_label has TAN-1..NAN-3, rename
sc_corr["neu_subtype"] = sc_corr["tan_label"].map(RENAME).fillna(sc_corr["tan_label"])
sc_corr = sc_corr[["neu_subtype", "MP", "spearman_rho", "p", "n", "padj_full"]].copy()
sc_corr = sc_corr.rename(columns={"p": "pvalue", "padj_full": "padj_BH"})
# drop rows for unrenamed labels (shouldn't exist but defensive)
sc_corr = sc_corr.dropna(subset=["neu_subtype"])
sc_corr["neu_subtype"] = pd.Categorical(sc_corr["neu_subtype"], categories=ORDER, ordered=True)
sc_corr = sc_corr.sort_values(["neu_subtype", "MP"]).reset_index(drop=True)
assert_clean(sc_corr, "5F")
sc_corr.to_csv(OUT / "fig5f_neu_mp_correlation.csv", index=False)
print(f"  {sc_corr.shape} → fig5f_neu_mp_correlation.csv")

# also the per-patient combined matrix (for scatter regeneration in R)
combined = pd.read_csv(RES / "step26_patient_combined.csv")
combined.columns = [c.replace("frac_TAN-1", "frac_Neu_Inflammatory")
                     .replace("frac_TAN-2", "frac_Neu_IFN_response")
                     .replace("frac_TAN-3", "frac_Neu_Angiogenic")
                     .replace("frac_TAN-4", "frac_Neu_Metastatic")
                     .replace("frac_NAN-1", "frac_Neu_ECM_remodeling")
                     .replace("frac_NAN-2", "frac_Neu_OSM_priming")
                     .replace("frac_NAN-3", "frac_Neu_OSM_low")
                     .replace("frac_Neutrophils", "frac_Neu_unclassified")
                    for c in combined.columns]
combined = combined.rename(columns={"Unnamed: 0": "sample_id"})
assert_clean(combined, "5F-scatter")
combined.to_csv(OUT / "fig5f_patient_scatter_data.csv", index=False)
print(f"  {combined.shape} → fig5f_patient_scatter_data.csv")

# 5G: LIANA top-30 LR pairs Neu_Inflammatory → Malignant
print("\n[5G] LIANA top-30")
# the renamed CSV from step25h
liana_inflam = pd.read_csv(RES / "step27b_liana_top30_Neu_Inflammatory_to_malignant.csv")
# columns: source, target, ligand_complex, receptor_complex, lr_means, magnitude_rank, etc.
# add pathway annotation
def lr_pathway(row):
    lig = str(row.get("ligand_complex", row.get("ligand", "")))
    rec = str(row.get("receptor_complex", row.get("receptor", "")))
    L = lig.upper()
    if "TGFB" in L: return "TGFb"
    if L in {"IL1B", "IL1A"} or "IL1" in L: return "IL1"
    if L == "OSM" or "OSM" in L: return "OSM"
    if L == "CXCL8" or "CXCL8" in L: return "CXCL8"
    if L.startswith("CXCL"): return "CXCL_other"
    if L.startswith("CCL"): return "CCL"
    if L == "PLAU" or "PLAU" in L: return "PLAU"
    if L == "VEGFA" or "VEGF" in L: return "VEGFA"
    if L == "TNF": return "TNF"
    if L in {"AREG", "EREG", "HBEGF", "TGFA", "EPGN", "EGF"}: return "EGFR_ligand"
    if L in {"SPP1"}: return "SPP1"
    if L in {"FN1"}: return "FN1"
    return "Other"

# canonical column names for liana 1.7
lig_col = "ligand_complex" if "ligand_complex" in liana_inflam.columns else "ligand"
rec_col = "receptor_complex" if "receptor_complex" in liana_inflam.columns else "receptor"
key_score = "magnitude_rank" if "magnitude_rank" in liana_inflam.columns else "lr_means"
key_spec = "specificity_rank" if "specificity_rank" in liana_inflam.columns else None

liana_inflam["pathway"] = liana_inflam.apply(lr_pathway, axis=1)
df_g = pd.DataFrame({
    "sender":            liana_inflam.get("source_renamed", liana_inflam["source"]).astype(str),
    "receiver":          liana_inflam.get("target_renamed", liana_inflam["target"]).astype(str),
    "ligand":            liana_inflam[lig_col].astype(str),
    "receptor":          liana_inflam[rec_col].astype(str),
    "magnitude_rank":    liana_inflam[key_score].astype(float),
    "specificity_rank":  liana_inflam[key_spec].astype(float) if key_spec and key_spec in liana_inflam.columns else np.nan,
    "lr_means":          liana_inflam.get("lr_means", pd.Series([np.nan]*len(liana_inflam))).astype(float),
    "expr_prop":         liana_inflam.get("expr_prop", pd.Series([np.nan]*len(liana_inflam))).astype(float),
    "pathway":           liana_inflam["pathway"].astype(str),
}).head(30)
# strip the 'Mal_' prefix from receiver for cleaner R labels
df_g["receiver_short"] = df_g["receiver"].str.replace("Mal_", "")
assert_clean(df_g, "5G")
df_g.to_csv(OUT / "fig5g_liana_lr_pairs.csv", index=False)
print(f"  {df_g.shape} → fig5g_liana_lr_pairs.csv (top-30 Neu_Inflammatory → Mal_*)")

# also build the combined focus heatmap data (multi-sender × Mal_MP*) for chord/heatmap
focus = pd.read_csv(RES / "step27b_liana_focus_pathways_renamed.csv")
focus["pathway_simple"] = focus.apply(lr_pathway, axis=1)
sender_col = "source_renamed" if "source_renamed" in focus.columns else "source"
focus_out = pd.DataFrame({
    "sender":           focus[sender_col].astype(str),
    "receiver":         focus.get("target_renamed", focus["target"]).astype(str),
    "ligand":           focus[lig_col].astype(str),
    "receptor":         focus[rec_col].astype(str),
    "magnitude_rank":   focus[key_score].astype(float),
    "pathway":          focus["pathway"].astype(str) if "pathway" in focus.columns else focus["pathway_simple"],
})
focus_out["receiver_short"] = focus_out["receiver"].str.replace("Mal_", "")
focus_out["sender_short"] = focus_out["sender"].str.replace("Neu_", "")
focus_out = focus_out.dropna(subset=["sender", "receiver"]).reset_index(drop=True)
assert_clean(focus_out, "5G-focus")
focus_out.to_csv(OUT / "fig5g_liana_focus_all_senders.csv", index=False)
print(f"  {focus_out.shape} → fig5g_liana_focus_all_senders.csv")

# 5H: TCGA correlation matrix (long format with rho + p)
print("\n[5H] TCGA correlation matrix")
tcga_corr = pd.read_csv(RES / "step28_tcga_correlation_matrix.csv")
tcga_corr = tcga_corr.rename(columns={"sig": "neu_subtype", "p": "pvalue"})
tcga_corr.to_csv(OUT / "fig5h_tcga_correlation_matrix.csv", index=False)
# wide form (subtype × MP) for heatmap
tcga_wide = tcga_corr.pivot_table(index="neu_subtype", columns="MP", values="spearman_rho")
tcga_p_wide = tcga_corr.pivot_table(index="neu_subtype", columns="MP", values="pvalue")
tcga_wide_combined = tcga_wide.copy()
tcga_wide_combined.columns = [f"{c}_rho" for c in tcga_wide_combined.columns]
for col in tcga_p_wide.columns:
    tcga_wide_combined[f"{col}_pvalue"] = tcga_p_wide[col]
tcga_wide_combined = tcga_wide_combined.reset_index()
assert_clean(tcga_wide_combined, "5H-wide")
tcga_wide_combined.to_csv(OUT / "fig5h_tcga_correlation_wide.csv", index=False)
assert_clean(tcga_corr, "5H")
print(f"  long {tcga_corr.shape} + wide {tcga_wide_combined.shape}")

# 5I: KM 4-group data (TCGA Neu_Metastatic × MP2)
print("\n[5I] KM 4-group data")
# Need clinical OS + per-sample Neu_Metastatic + MP2 score
# Reload TCGA combined scores + clinical
combined_scores = pd.read_csv(RES / "step28_tcga_combined_scores.csv.gz", index_col=0)
clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin_pt = clin[clin["sample_type"] == "Primary Tumor"].copy()
df_i = clin_pt.set_index("sample_barcode").join(combined_scores, how="inner").reset_index()
df_i["event"] = (df_i["vital_status"].str.strip().str.lower() == "dead").astype(int)
df_i["time_days"] = np.where(df_i["event"] == 1, df_i["days_to_death"], df_i["days_to_last_follow_up"])
df_i = df_i[df_i["time_days"].notna() & (df_i["time_days"] > 0)].copy()

# 4-group classification
sm = df_i["Neu_Metastatic"].median()
mm = df_i["MP2"].median()
df_i["neu_met_group"] = np.where(df_i["Neu_Metastatic"] >= sm, "High", "Low")
df_i["mp2_group"]     = np.where(df_i["MP2"] >= mm, "High", "Low")
df_i["combined_group"] = df_i["neu_met_group"] + "_NeuMet/" + df_i["mp2_group"] + "_MP2"
# clean ordered factor
combo_levels = ["Low_NeuMet/Low_MP2", "Low_NeuMet/High_MP2",
                "High_NeuMet/Low_MP2", "High_NeuMet/High_MP2"]
df_i["combined_group"] = pd.Categorical(df_i["combined_group"], categories=combo_levels, ordered=True)

# select clean columns
out_cols_5i = ["sample_barcode", "patient_id" if "patient_id" in df_i.columns else None,
               "neu_met_group", "mp2_group", "combined_group",
               "Neu_Metastatic", "MP2", "MP1", "MP3", "MP4",
               "time_days", "event", "vital_status", "ajcc_stage", "gender",
               "age_at_diagnosis"]
out_cols_5i = [c for c in out_cols_5i if c is not None and c in df_i.columns]
df_i_out = df_i[out_cols_5i].copy()
df_i_out = df_i_out.rename(columns={"sample_barcode": "sample_id",
                                     "time_days": "OS_days",
                                     "event": "OS_status",
                                     "ajcc_stage": "stage"})
assert_clean(df_i_out, "5I")
df_i_out.to_csv(OUT / "fig5i_km_data.csv", index=False)
print(f"  {df_i_out.shape} → fig5i_km_data.csv (n events={df_i_out['OS_status'].sum()})")

# also build the parallel 4-group tables for the other axes (supp)
for sig, mp, suffix in [("Neu_Inflammatory", "MP1", "neuInflam_mp1"),
                         ("Neu_OSM_priming", "MP1", "neuOSM_mp1"),
                         ("Neu_ECM_remodeling", "MP1", "neuECM_mp1")]:
    s_med = df_i[sig].median(); m_med = df_i[mp].median()
    g1 = np.where(df_i[sig] >= s_med, "High", "Low")
    g2 = np.where(df_i[mp] >= m_med, "High", "Low")
    combo = [f"{a}_{sig.replace('Neu_','')}/{b}_{mp}" for a, b in zip(g1, g2)]
    sub = df_i[["sample_barcode", sig, mp, "time_days", "event"]].copy()
    sub.columns = ["sample_id", sig, mp, "OS_days", "OS_status"]
    sub["sig_group"] = g1
    sub["mp_group"] = g2
    sub["combined_group"] = combo
    sub.to_csv(OUT / f"fig5i_km_{suffix}.csv", index=False)
    print(f"  + fig5i_km_{suffix}.csv ({sub.shape})")

# 5J: Cox forest data (univariate + multivariate, combined table)
print("\n[5J] Cox forest data")
uni = pd.read_csv(RES / "step28_cox_univariate.csv")
uni["model"] = "univariate"
uni["covariate"] = ""
uni = uni.rename(columns={"signature": "variable",
                          "HR_lo": "CI_lower", "HR_hi": "CI_upper",
                          "p": "pvalue"})
multi = pd.read_csv(RES / "step28_cox_multivariate.csv")
multi = multi.rename(columns={"HR_lo": "CI_lower", "HR_hi": "CI_upper", "p": "pvalue"})
# columns from step28: variable, coef, HR, CI_lower, CI_upper, pvalue, model, n
multi["covariate"] = multi["model"]
multi["model"] = "multivariate"
# unify schema
common = ["variable", "HR", "CI_lower", "CI_upper", "pvalue", "n", "model", "covariate"]
for df in (uni, multi):
    for c in common:
        if c not in df.columns:
            df[c] = np.nan
df_j = pd.concat([uni[common], multi[common]], ignore_index=True)

# group flag for HR direction
df_j["effect"] = np.where((df_j["HR"] > 1) & (df_j["pvalue"] < 0.05), "risk",
                  np.where((df_j["HR"] < 1) & (df_j["pvalue"] < 0.05), "protective", "n.s."))

assert_clean(df_j, "5J")
df_j.to_csv(OUT / "fig5j_cox_forest_data.csv", index=False)
print(f"  {df_j.shape} → fig5j_cox_forest_data.csv (uni n={(df_j['model']=='univariate').sum()}, multi n={(df_j['model']=='multivariate').sum()})")

# Bonus: panel-color reference (so R can match)
print("\n[bonus] color palette + manifest")
PALETTE = {
    "Neu_Inflammatory":   "#d62728",
    "Neu_Angiogenic":     "#ff7f0e",
    "Neu_Metastatic":     "#9467bd",
    "Neu_IFN_response":   "#bcbd22",
    "Neu_OSM_priming":    "#1f77b4",
    "Neu_OSM_low":        "#17becf",
    "Neu_ECM_remodeling": "#2ca02c",
    "Neu_unclassified":   "#7f7f7f",
}
pd.DataFrame([{"neu_subtype": k, "color_hex": v, "display_order": i+1}
              for i, (k, v) in enumerate(PALETTE.items())]).to_csv(
    OUT / "fig5_palette.csv", index=False)

# Manifest
manifest = []
for f in sorted(OUT.glob("*.csv*")):
    manifest.append({"file": f.name, "size_kb": round(f.stat().st_size / 1024, 1)})
pd.DataFrame(manifest).to_csv(OUT / "_manifest.csv", index=False)

print(f"\nelapsed: {time.time()-t0:.1f}s")
print("DONE.")
