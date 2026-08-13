"""Update Supplementary Tables S1-S11 with Q/A audit fixes."""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

OUT = Path("${WORK_ROOT}/Supplementary_Tables")


def write_titled(xw, sheet, title, subtitle, df):
    rows = [[title] + [None] * max(0, len(df.columns) - 1)]
    if subtitle:
        rows.append([subtitle] + [None] * max(0, len(df.columns) - 1))
    rows.append([None] * len(df.columns))
    rows.append(list(df.columns))
    rows.extend(df.values.tolist())
    pd.DataFrame(rows).to_excel(xw, sheet_name=sheet[:31], index=False, header=False)


def read(path, **kw):
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_csv(p, **kw)


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def stars(p):
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"


# S1 - drop GSE127465
print("[S1] data sources")
data_sources = pd.DataFrame([
    ("scRNA-seq",  "GSE131907", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE143423", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143423", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE148071", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE164789", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164789", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE189357", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189357", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE253013", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE123902", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902", "Integration analysis", "10x"),
    ("scRNA-seq",  "Salcher 2022 LUAD atlas", "https://doi.org/10.1016/j.ccell.2022.10.008", "Neutrophil / myeloid high-resolution reference", "10x / mixed"),
    ("Spatial RNA-seq", "E-MTAB-13530", "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530", "Spatial validation (12 sections)", "10x Visium"),
    ("Spatial RNA-seq", "Takano 2024 (Visium FF)", "https://doi.org/10.1038/s41467-024-54671-7", "Spatial validation (8 fresh-frozen sections)", "10x Visium"),
    ("Bulk RNA-seq",   "TCGA-LUAD", "https://portal.gdc.cancer.gov/projects/TCGA-LUAD", "Survival, Cox, KM (Fig 3, Fig 8O/P)", "Illumina"),
    ("Bulk RNA-seq",   "GSE207422", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422", "Anti-PD-(L)1 neoadjuvant response (Fig 8M)", "Illumina"),
    ("Bulk RNA-seq",   "GSE126044", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044", "Anti-PD-1 R vs NR (Fig 8N)", "Illumina"),
    ("scRNA-seq (validation)", "GSE135222", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222", "Single-cell anti-PD-1 R vs NR validation (Fig S11D-F)", "10x"),
    ("CRISPR essentiality", "DepMap 24Q2", "https://depmap.org/portal/data_page/?tab=allData", "LUAD cell-line lethality benchmark (Fig 8B/G)", "Avana CRISPR"),
], columns=["DataType", "DataSource", "Accession Link", "Application purpose", "Platform"])

# Sample info from cell metadata
meta = read("${WORK_ROOT}/luad_figures/fig_s1/s1_cell_metadata.csv.gz",
            usecols=["dataset","sample_id","patient_id","tissue_type"])
sample_info = (meta.drop_duplicates(["dataset","sample_id"])
                   .rename(columns={"dataset":"Dataset","sample_id":"Sample",
                                    "patient_id":"Patient","tissue_type":"Tissue"})
                   .sort_values(["Dataset","Sample"])
                   .reset_index(drop=True))
spatial_samples = pd.DataFrame([
    ("E-MTAB-13530", s, p, "Tumor / boundary") for s, p in [
        ("P10_T1","P10"),("P10_T2","P10"),("P10_T3","P10"),("P10_T4","P10"),
        ("P15_T1","P15"),("P15_T2","P15"),("P16_T1","P16"),("P16_T2","P16"),
        ("P24_T1","P24"),("P24_T2","P24"),("P25_T1","P25"),("P25_T2","P25"),
    ]
] + [
    ("Takano 2024 (FF)", f"LUAD_No_{n}", f"Patient_{n}", "Invasive LUAD")
    for n in (1, 2, 3, 4, 5, 14, 16, 17)
], columns=["Dataset","Sample","Patient","Tissue"])
sample_info_full = pd.concat([sample_info, spatial_samples], ignore_index=True)

with pd.ExcelWriter(OUT / "Supplementary_Table_S1_Data_Sources.xlsx") as xw:
    write_titled(xw, "Data sources",
        "Table S1. Data sources used in this study.",
        "All cohorts are public; GSE127465 was screened during pilot work and is not used in the final analyses.",
        data_sources)
    write_titled(xw, "Sample info",
        "Table S1 (continued). Per-sample metadata for the 256 single-cell samples and 20 spatial sections.",
        None,
        sample_info_full)


# S2 - drop empty clusters + cluster->MP mapping
print("[S2] cNMF top genes")
gep_top = read("${WORK_ROOT}/luad_figures/fig_s2/s2_gep_top30_genes.csv")
corr = read("${WORK_ROOT}/luad_figures/fig_s2/gep_spearman_corr.csv", index_col=0)
mp_map = read("${WORK_ROOT}/luad_figures/fig2/cnmf_consensus_mp_annotation.csv")

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
keep = [c for c in corr.index if c in corr.columns]
sub = corr.loc[keep, keep].values
sub = (sub + sub.T) / 2
np.fill_diagonal(sub, 1.0)
dist = np.clip(1 - sub, 0, None)
np.fill_diagonal(dist, 0.0)
Z = linkage(squareform(dist, checks=False), method="average")
clust15 = fcluster(Z, t=15, criterion="maxclust")
gep_to_clust = pd.Series(clust15, index=keep)

# Compute cluster -> MP majority assignment
# `cnmf_consensus_mp_annotation.csv` already maps program_id (cNMF_X) -> MP
cluster_to_mp = {}
if mp_map is not None and {"program_id", "MP"}.issubset(mp_map.columns):
    for _, r in mp_map.iterrows():
        m = re.match(r"cNMF_(\d+)", str(r["program_id"]))
        if m:
            cluster_to_mp[int(m.group(1))] = r["MP"]

gep_top_in = gep_top[gep_top["gep_id"].isin(gep_to_clust.index)].copy()
gep_top_in["cluster"] = gep_to_clust.loc[gep_top_in["gep_id"]].values
TOP_N = 30
out_cols = {}
header_map = {}
for cid in range(1, int(clust15.max()) + 1):
    sub_df = gep_top_in[gep_top_in["cluster"] == cid]
    if len(sub_df) == 0:
        continue   # drop empty clusters
    top = (sub_df.groupby("gene")["rank"].mean()
                  .sort_values()
                  .head(TOP_N).index.tolist())
    while len(top) < TOP_N:
        top.append(None)
    label = f"cNMF_{cid}"
    if cid in cluster_to_mp:
        label = f"cNMF_{cid} ({cluster_to_mp[cid]})"
    out_cols[label] = top
    header_map[cid] = label
s2 = pd.DataFrame(out_cols)

mp_caption = ("Table S2. Top-30 marker genes per cNMF cluster.")
mp_subtitle = (
    "k = 15 hierarchical cuts of the 770 program-level GEP correlation matrix; "
    "empty clusters (insufficient stable members) are omitted. Column header "
    "shows cluster index and majority MP assignment in parentheses."
)
with pd.ExcelWriter(OUT / "Supplementary_Table_S2_cNMF_Top_Genes.xlsx") as xw:
    write_titled(xw, "cNMF_top_genes", mp_caption, mp_subtitle, s2)


# S3 - MP DEGs caption update (MP1-MP5)
print("[S3] MP DEGs")
mp_sig = read("${WORK_ROOT}/luad_figures/fig2/mp_signatures_top100.csv")
mp_sig_out = (mp_sig.rename(columns={
        "MP": "cluster",
        "rank": "rank_within_MP",
        "consensus_score": "consensus_score",
        "n_GEPs_in_MP": "n_GEPs_supporting",
    })[["cluster","rank_within_MP","gene","consensus_score","n_GEPs_supporting"]])
with pd.ExcelWriter(OUT / "Supplementary_Table_S3_Meta_Program_DEGs.xlsx") as xw:
    write_titled(xw, "Meta Programs DEGs",
        "Table S3. Meta-Program (MP1-MP5) signature genes - top 100 per MP.",
        ("MP1-MP4 are the validated programs reported in the main text; MP5 is "
         "included for completeness but represents single-patient noise and is "
         "excluded from downstream analyses. consensus_score quantifies how "
         "consistently a gene is selected as a top-30 program member across "
         "patient GEPs assigned to that MP."),
        mp_sig_out)


# S4 - Neutrophil subtype markers (no caption change)
print("[S4] Neutrophil subtype DEGs")
neu = read("${WORK_ROOT}/luad_figures/fig5/data/fig5c_canonical_markers.csv")
neu_out = (neu.rename(columns={
        "neu_subtype": "cluster",
        "mean_expression": "mean_log_expr",
    })[["cluster","gene","mean_log_expr","pct_expressing","canonical_for"]])
with pd.ExcelWriter(OUT / "Supplementary_Table_S4_Neutrophil_Subtype_DEGs.xlsx") as xw:
    write_titled(xw, "Neutrophil subtype DEGs",
        "Table S4. Neutrophil subtype canonical marker expression.",
        ("Source: Salcher 2022 LUAD atlas + integrated scRNA-seq (Figure 5C). "
         "mean_log_expr = mean log1p-normalised expression in the subtype; "
         "pct_expressing = percentage of cells in the subtype with expr > 0; "
         "canonical_for = subtype the gene was originally curated as a marker of."),
        neu_out)


# S5 - rename sheet + correct scope (myeloid not just macrophage)
print("[S5] Myeloid subtype DEGs")
m1 = read("${WORK_ROOT}/luad_figures/fig4/panel_GM_subset_markers.csv")
m1_out = (m1.rename(columns={
        "subtype": "cluster",
        "mean_log1p": "mean_log_expr",
        "panel_origin": "canonical_for",
    })[["cluster","gene","mean_log_expr","pct_expressing","canonical_for"]])
m2 = read("${WORK_ROOT}/luad_figures/fig4/panel_N_spp1_vs_c1qc_deg.csv")
with pd.ExcelWriter(OUT / "Supplementary_Table_S5_Myeloid_Subtype_DEGs.xlsx") as xw:
    write_titled(xw, "Myeloid subtype markers",
        "Table S5a. Myeloid subtype marker expression - all 13 myeloid subtypes (panel_GM).",
        ("Source: integrated scRNA-seq myeloid compartment (Figure 4). Includes "
         "Macro_*, Mono_nonclassical, Neutrophil and dendritic-cell subtypes. "
         "mean_log_expr = mean log1p expression; canonical_for = subset that the "
         "gene was originally curated for."),
        m1_out)
    if m2 is not None:
        m2_out = m2.copy()
        m2_out.columns = ["gene","score","avg_log2FC","p_val","p_val_adj"][:len(m2_out.columns)]
        write_titled(xw, "SPP1 vs C1QC pseudobulk DEG",
            "Table S5b. Macro_SPP1 vs Macro_C1QC pseudobulk DEGs.",
            "Pseudobulk DESeq2 (one row per patient x subtype, Figure 4N).",
            m2_out)


# S6 - caption clarifying Sig=1 vs final pool filter
print("[S6] Geneformer KO transitions")
GF = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/geneformer_500_stats")
sheet_specs = [
    ("Macro-SPP1-to-Macro-C1QC", GF / "macro_spp1_to_c1qc_stats.csv"),
    ("Mal-MP3-to-MP1",           GF / "mal_mp3_to_mp1_stats.csv"),
    ("Neu-OSM-priming-to-low",   GF / "neu_osm_priming_to_low_stats.csv"),
]
with pd.ExcelWriter(OUT / "Supplementary_Table_S6_Geneformer_KO_Transitions.xlsx") as xw:
    for sn, p in sheet_specs:
        if not p.exists():
            continue
        df = pd.read_csv(p)
        if df.columns[0].startswith("Unnamed") or df.columns[0] == "":
            df = df.drop(columns=df.columns[0])
        df = df[["Gene","Gene_name","Ensembl_ID","Shift_to_goal_end",
                 "Goal_end_vs_random_pval","Goal_end_FDR","N_Detections","Sig"]]
        df = df.sort_values("Shift_to_goal_end", ascending=False)
        write_titled(xw, sn,
            f"Table S6. Geneformer KO transition: {sn}.",
            ("Geneformer V2 (104 M parameters) in-silico KO/over-expression "
             "simulation (Figure 8 perturbation analysis). Shift_to_goal_end = "
             "mean shift toward the desired endpoint embedding (positive = "
             "closer to goal). Sig = 1 indicates Goal_end_FDR < 0.05 against a "
             "random-gene null. The Fig 8A candidate pool is further filtered "
             "to Shift_to_goal_end > 0 AND N_Detections >= 5 to remove low-"
             "confidence calls."),
            df)


# S7 - caption corrections + risk-score Cox addition
print("[S7] Survival statistics")
fig3 = Path("${WORK_ROOT}/luad_figures/fig3")
fig8d = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data")

cox_uv = read(fig3 / "tcga_luad_mp_cox_univariate.csv")
cox_mv = read(fig3 / "tcga_luad_mp_cox_multivariate.csv")
km_med = read(fig3 / "tcga_luad_mp_km_logrank.csv")
km_opt = read(fig3 / "tcga_optimal_cutoff_km.csv")
km_risk = read(fig3 / "tcga_risk_score_km.csv")
km_srsf9 = read(fig8d / "8O_km_SRSF9_stats.csv")
km_sec61g = read(fig8d / "8P_km_SEC61G_stats.csv")
risk_cox = read(fig3 / "tcga_risk_score_cox.csv")

with pd.ExcelWriter(OUT / "Supplementary_Table_S7_Survival_Statistics.xlsx") as xw:
    if cox_uv is not None:
        write_titled(xw, "Cox univariate",
            "Table S7a. TCGA-LUAD Cox univariate regression for each MP score (Fig 3A forest plot input).",
            "MP1-4 z-scored across the cohort; outcome = overall survival.",
            cox_uv)
    if cox_mv is not None:
        write_titled(xw, "Cox multivariate",
            "Table S7b. TCGA-LUAD Cox multivariate regression (Age + Stage + MP1-4) - Fig 3A forest plot.",
            "Covariates Age (continuous), Stage (ordinal I-IV), MP1-MP4 z-scored.",
            cox_mv)
    if km_med is not None:
        write_titled(xw, "KM median split",
            "Table S7c. TCGA-LUAD KM log-rank by median MP score split (Fig 3B).",
            None, km_med)
    if km_opt is not None:
        write_titled(xw, "KM optimal cutoff",
            "Table S7d. TCGA-LUAD KM with optimal-cutoff (max-rank) split per MP (Fig 3D).",
            "Optimal cutoff selected to maximise log-rank statistic (constrained to 30-70 percentile).",
            km_opt)
    if km_risk is not None:
        write_titled(xw, "Risk score KM",
            "Table S7e. TCGA-LUAD risk-score KM (Fig 3F).",
            "Risk score = z(MP2) - z(MP4) per multivariate Cox; tertile split.",
            km_risk)
    pieces = []
    for stat, gene in [(km_sec61g, "SEC61G"), (km_srsf9, "SRSF9")]:
        if stat is None:
            continue
        s = stat.copy()
        if "gene" in s.columns:
            s = s.drop(columns=["gene"])
        s.insert(0, "gene", gene)
        pieces.append(s)
    if pieces:
        lg = pd.concat(pieces, ignore_index=True)
        write_titled(xw, "Lead gene KM",
            "Table S7f. TCGA-LUAD KM for the two lead genes SEC61G and SRSF9 (Fig 8O/P).",
            ("ANGPTL4 was assessed but did not reach significance "
             "(log-rank p > 0.05) and is not shown."),
            lg)
    if risk_cox is not None:
        write_titled(xw, "Risk score Cox",
            "Table S7g. Cox univariate and multivariate hazard ratios for the composite risk score (Fig 3F).",
            "Risk score treated as a continuous covariate; multivariate model adjusts for Age and Stage.",
            risk_cox)


# S8 - cohort label fix + Takano lead-gene per-section computation
print("[S8] Spatial validation")
r_data = Path("${DATA_ROOT}/ST/results/r_data")
roi_agg = read(r_data / "roi_vs_nonroi_aggregate_pvalues.csv")
roi_sec = read(r_data / "roi_vs_nonroi_stats_with_pvalues.csv")
roi_consist = read(r_data / "roi_per_section_consistency.csv")
misty = read(r_data / "misty_aggregated_importance.csv")
sec_mp_e = read(fig8d / "8I_section_MP_scores.csv")
sec_mp_o = read(fig8d / "okamura_section_MP_scores.csv")  # may exist
lead_sec_e = read(fig8d / "per_section_lead_gene_sig.csv")
lead_coh_e = read(fig8d / "per_cohort_lead_gene_sig.csv")

# rename Okamura -> Takano 2024 in any frame that has a cohort column
for df in (roi_agg, roi_sec, lead_sec_e, lead_coh_e, misty):
    if df is not None and "cohort" in df.columns:
        df["cohort"] = df["cohort"].replace({
            "Okamura": "Takano 2024",
            "EMTAB13530": "E-MTAB-13530",
        })
# MISTy table sometimes uses a 'dataset' column
if misty is not None and "dataset" in misty.columns:
    misty["dataset"] = misty["dataset"].replace({
        "Okamura": "Takano 2024",
        "EMTAB13530": "E-MTAB-13530",
    })

# ---- Takano lead-gene computation ----
print("  computing Takano lead gene Mann-Whitney from cohort.h5ad + per_section CSVs")
TAKANO_FF = [f"LUAD_No_{n}" for n in (1, 2, 3, 4, 5, 14, 16, 17)]
GENES = ["SEC61G", "SRSF9", "ANGPTL4"]
takano_h5 = "${DATA_ROOT}/ST/results/step09_okamura_validation/cohort.h5ad"
ad = sc.read_h5ad(takano_h5)
keep_mask = ad.obs["sample"].isin(TAKANO_FF)
ad = ad[keep_mask].copy()
# build per-section data merging MP3 + ct_Malignant from per_section CSV
takano_per_sec_rows = []
takano_buf_for_cohort = {g: {"a": [], "b": []} for g in GENES}
for sample in TAKANO_FF:
    csv_p = f"${DATA_ROOT}/ST/results/r_data/per_section/Okamura__{sample}.csv"
    if not Path(csv_p).exists():
        continue
    sec_df = pd.read_csv(csv_p)
    sec_idx = sec_df.set_index("spot_id")
    obs_in_section = ad[ad.obs["sample"] == sample].copy()
    common = obs_in_section.obs_names.intersection(sec_idx.index)
    if len(common) < 20:
        continue
    obs_sub = obs_in_section[common].copy()
    sec_aligned = sec_idx.loc[common]
    mp3 = sec_aligned["MP3_score"].values.astype(float)
    mal = sec_aligned["ct_Malignant"].values.astype(float)
    z_mp3 = (mp3 - mp3.mean()) / (mp3.std(ddof=0) + 1e-12)
    z_mal = (mal - mal.mean()) / (mal.std(ddof=0) + 1e-12)
    new_roi = (z_mp3 > 0.5) & (z_mal > 0.5)
    if new_roi.sum() < 3 or (~new_roi).sum() < 3:
        continue
    sec_buf = []
    for g in GENES:
        if g not in obs_sub.var_names:
            continue
        e = obs_sub[:, g].X
        e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
        a = e[new_roi]; b = e[~new_roi]
        if len(a) < 3 or len(b) < 3:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        sec_buf.append({
            "cohort": "Takano 2024",
            "sample": sample,
            "gene": g,
            "n_roi": int(len(a)),
            "n_nonroi": int(len(b)),
            "mean_roi": float(a.mean()),
            "mean_nonroi": float(b.mean()),
            "delta": float(a.mean() - b.mean()),
            "U_stat": float(u),
            "p_raw": float(p),
        })
        takano_buf_for_cohort[g]["a"].extend(a.tolist())
        takano_buf_for_cohort[g]["b"].extend(b.tolist())
    if sec_buf:
        q = bh_fdr(np.asarray([r["p_raw"] for r in sec_buf]))
        for row, qi in zip(sec_buf, q):
            row["p_fdr"] = float(qi)
            row["sig"]   = stars(qi)
            takano_per_sec_rows.append(row)

# Cohort-level test for Takano
takano_cohort_rows = []
for g in GENES:
    a = np.asarray(takano_buf_for_cohort[g]["a"])
    b = np.asarray(takano_buf_for_cohort[g]["b"])
    if len(a) < 5 or len(b) < 5:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    takano_cohort_rows.append({
        "cohort": "Takano 2024",
        "gene": g,
        "n_roi": int(len(a)),
        "n_nonroi": int(len(b)),
        "mean_roi": float(a.mean()),
        "mean_nonroi": float(b.mean()),
        "delta": float(a.mean() - b.mean()),
        "U_stat": float(u),
        "p_raw": float(p),
    })
if takano_cohort_rows:
    q = bh_fdr(np.asarray([r["p_raw"] for r in takano_cohort_rows]))
    for row, qi in zip(takano_cohort_rows, q):
        row["p_fdr"] = float(qi); row["sig"] = stars(qi)

# concatenate with EMTAB rows from earlier
lead_sec_combined = pd.concat([lead_sec_e, pd.DataFrame(takano_per_sec_rows)], ignore_index=True) if lead_sec_e is not None else pd.DataFrame(takano_per_sec_rows)
lead_coh_combined = pd.concat([lead_coh_e, pd.DataFrame(takano_cohort_rows)], ignore_index=True) if lead_coh_e is not None else pd.DataFrame(takano_cohort_rows)

# also harmonise the cohort label in the EMTAB part of the combined frames
for df in (lead_sec_combined, lead_coh_combined):
    if "cohort" in df.columns:
        df["cohort"] = df["cohort"].replace({
            "Okamura": "Takano 2024",
            "EMTAB13530": "E-MTAB-13530",
        })

# Also build an "S8h Per-section MP scores - Takano" frame
takano_mp_rows = []
for sample in TAKANO_FF:
    csv_p = f"${DATA_ROOT}/ST/results/r_data/per_section/Okamura__{sample}.csv"
    if not Path(csv_p).exists():
        continue
    sec_df = pd.read_csv(csv_p)
    cols = [c for c in ["MP1_score","MP2_score","MP3_score","MP4_score",
                         "MP5_score","ct_Malignant","ct_Neu_Inflammatory",
                         "ct_Macro_SPP1"] if c in sec_df.columns]
    means = sec_df[cols].mean()
    row = {"cohort": "Takano 2024", "sample": sample}
    row.update({k: float(v) for k, v in means.items()})
    takano_mp_rows.append(row)
takano_mp = pd.DataFrame(takano_mp_rows)

with pd.ExcelWriter(OUT / "Supplementary_Table_S8_Spatial_Validation_Statistics.xlsx") as xw:
    if roi_agg is not None:
        write_titled(xw, "ROI vs non-ROI cohort",
            "Table S8a. ROI vs non-ROI per-cohort statistics (Fig 7Q, Fig S10J).",
            ("Per-spot Mann-Whitney U test (two-sided) pooled within cohort across "
             "sections; BH-FDR within cohort across all 67 metrics."),
            roi_agg)
    if roi_sec is not None:
        write_titled(xw, "ROI vs non-ROI per-section",
            "Table S8b. ROI vs non-ROI per-section statistics (Fig S10K input).",
            ("Mann-Whitney U test per (sample, metric); BH-FDR within each section "
             "across metrics. Sections with < 5 ROI or < 5 non-ROI spots dropped."),
            roi_sec)
    if roi_consist is not None:
        write_titled(xw, "Per-section consistency",
            "Table S8c. Per-metric consistency summary across all sections.",
            ("n_sig_strong = sections with FDR<0.01; n_sig_any = sections with "
             "FDR<0.05; sign_consistent = sections agreeing with the metric's "
             "median delta sign."),
            roi_consist)
    if misty is not None:
        write_titled(xw, "MISTy importance",
            "Table S8d. MISTy cell-type interaction importances (Fig 7H).",
            ("Aggregated importance per (target, predictor) across views (intra/"
             "juxta/para); higher = stronger spatial dependency. Cohort labels "
             "use Takano 2024 throughout."),
            misty)
    if sec_mp_e is not None:
        e = sec_mp_e.copy()
        if "cohort" in e.columns:
            e["cohort"] = e["cohort"].replace({"Okamura": "Takano 2024",
                                                "EMTAB13530": "E-MTAB-13530"})
        write_titled(xw, "Per-section MP scores E-MTAB",
            "Table S8e. Per-section mean MP1-4 scores - E-MTAB-13530 (Fig 8I auxiliary).",
            None, e)
    if not takano_mp.empty:
        write_titled(xw, "Per-section MP scores Takano",
            "Table S8h. Per-section mean MP1-5 / cell-type abundances - Takano 2024 cohort.",
            ("Computed from the Fig 7 per-section CSVs (8 fresh-frozen sections). "
             "MP3 / ct_Malignant define the tumor-intrinsic ROI used in Fig S11Q/R."),
            takano_mp)
    if not lead_coh_combined.empty:
        write_titled(xw, "Lead gene cohort sig",
            "Table S8f. Tumor-intrinsic ROI vs non-ROI per-cohort Mann-Whitney + FDR for the 3 lead genes (Fig S11Q).",
            ("ROI = z(Malignant) > 0.5 AND z(MP3) > 0.5; per-spot test pooled "
             "within cohort; BH-FDR within cohort across 3 genes. Both E-MTAB-"
             "13530 and Takano 2024 cohorts now reported."),
            lead_coh_combined)
    if not lead_sec_combined.empty:
        write_titled(xw, "Lead gene per-section sig",
            "Table S8g. Tumor-intrinsic ROI vs non-ROI per-section Mann-Whitney + FDR for the 3 lead genes (Fig S11R).",
            ("Per-section spot-level test; BH-FDR within each section across 3 "
             "genes. E-MTAB-13530 + Takano 2024 sections combined."),
            lead_sec_combined)


# S9 - unchanged content (just rebuild for consistency)
print("[S9] Treatment / clinical validation")
depmap = read(fig8d / "8B_depmap_stats.csv")
expr_eff = read(fig8d / "8G_expr_vs_effect.csv")
tcga_tvn = read(fig8d / "8D_tcga_TvN_stats.csv")
mpr = read(fig8d / "8M_GSE207422_stats.csv")
volcano = read(fig8d / "8N_GSE126044_volcano.csv")
gse135 = read(fig8d / "S11D_GSE135222_stats.csv")
with pd.ExcelWriter(OUT / "Supplementary_Table_S9_Treatment_Clinical_Validation.xlsx") as xw:
    if depmap is not None:
        write_titled(xw, "DepMap LUAD essentiality",
            "Table S9a. DepMap 24Q2 CRISPR essentiality of lead genes in LUAD cell lines (Fig 8B).",
            "Gene effect (Chronos score, more negative = more essential); LUAD cohort vs non-LUAD; one-sided Mann-Whitney.",
            depmap)
    if expr_eff is not None:
        write_titled(xw, "Expr x essentiality",
            "Table S9b. Per-cell-line gene expression x CRISPR effect for lead genes (Fig 8G scatter).",
            "log2(TPM+1) vs gene_effect; LUAD lines only; Spearman rho per gene.",
            expr_eff)
    if tcga_tvn is not None:
        write_titled(xw, "TCGA Tumor vs Normal",
            "Table S9c. TCGA-LUAD Tumor vs Normal differential expression for lead genes (Fig 8D-F).",
            "Wilcoxon two-sided per gene; log2FC = log2(mean Tumor TPM+1) - log2(mean Normal TPM+1).",
            tcga_tvn)
    if mpr is not None:
        write_titled(xw, "GSE207422 MPR vs NMPR",
            "Table S9d. GSE207422 (neoadjuvant chemo-IO) MPR vs NMPR per-gene statistics (Fig 8M).",
            "Per-gene Wilcoxon two-sided + AUC of MPR vs NMPR responder classification.",
            mpr)
    if volcano is not None:
        write_titled(xw, "GSE126044 R vs NR",
            "Table S9e. GSE126044 (anti-PD-1) Responder vs Non-Responder full DE (Fig 8N volcano).",
            ("Wilcoxon per-gene; log2FC = R - NR; -log10(p) on y-axis. Lead genes "
             "flagged via is_lead column. (Note: a small number of gene symbols "
             "are auto-converted to dates by Excel - underlying data is not "
             "altered.)"),
            volcano)
    if gse135 is not None:
        write_titled(xw, "GSE135222 single-cell",
            "Table S9f. GSE135222 single-cell anti-PD-1 R vs NR validation of lead genes (Fig S11D-F).",
            "Per-gene Wilcoxon on log-normalised expression in tumor compartment; R vs NR.",
            gse135)


# S10 - Hallmark + Gavish caption with MP5 footnote
print("[S10] Pathway / hallmark enrichment")
fig2 = Path("${WORK_ROOT}/luad_figures/fig2")
hall_nes = read(fig2 / "hallmark_nes_heatmap.csv")
hall_fdr = read(fig2 / "hallmark_fdr_heatmap.csv")
gavish = read(fig2 / "gavish_top_matches.csv")
gavish_overlap = read(fig2 / "gavish_overlap.csv")
prog_e = read("${DATA_ROOT}/ST/results/step05_progeny/per_sample_mean_progeny.csv")
prog_o = read("${DATA_ROOT}/ST/results/step09_okamura_validation/per_sample_mean_progeny.csv")
prog_pieces = []
if prog_e is not None:
    e2 = prog_e.copy(); e2.insert(0, "cohort", "E-MTAB-13530"); prog_pieces.append(e2)
if prog_o is not None:
    o2 = prog_o.copy(); o2.insert(0, "cohort", "Takano 2024");  prog_pieces.append(o2)
prog_all = pd.concat(prog_pieces, ignore_index=True) if prog_pieces else None

with pd.ExcelWriter(OUT / "Supplementary_Table_S10_Pathway_Enrichment.xlsx") as xw:
    if hall_nes is not None:
        write_titled(xw, "Hallmark NES",
            "Table S10a. Hallmark gene-set NES per Meta-Program (Fig 2F).",
            "GSEA NES on MP scores vs all-other; positive = enriched in MP. 50 Hallmark sets x MP1-4.",
            hall_nes)
    if hall_fdr is not None:
        write_titled(xw, "Hallmark FDR",
            "Table S10b. BH-adjusted p-values matching the NES table (Fig 2F stars overlay).",
            "FDR < 0.05 / 0.01 / 0.001 = * / ** / *** in the figure.",
            hall_fdr)
    if gavish is not None:
        write_titled(xw, "Gavish pan-cancer alignment",
            "Table S10c. Top Gavish 41-MP pan-cancer match per cNMF cluster (Fig 2I).",
            ("Cosine similarity between cNMF cluster mean profile and each "
             "Gavish 2023 pan-cancer MP signature. MP5 is included for "
             "completeness; it represents single-patient noise and is "
             "excluded from the main analyses."),
            gavish)
    if gavish_overlap is not None:
        write_titled(xw, "Gavish overlap",
            "Table S10d. Gene-set overlap matrix between cNMF clusters and Gavish 2023 MPs (Fig 2I support).",
            "Overlap = |intersect| / min(|set_A|, |set_B|) using top-30 genes of each program.",
            gavish_overlap)
    if prog_all is not None:
        write_titled(xw, "PROGENy per-section",
            "Table S10e. Mean PROGENy 14-pathway activity per spatial section (Fig 7I, Fig S10H).",
            "Mean MLM-derived score per section across 14 pathways; both ST cohorts pooled.",
            prog_all)


# S11 - LIANA caption + COMMOT cohort label
print("[S11] TF / LR / COMMOT")
tf_top = read(fig3 / "tf_state_specific_top15.csv")
tf_z = read(fig3 / "tf_activity_mp_zscore.csv")
gs = read(fig3 / "geneswitches_results.csv")
gs_top = read(fig3 / "geneswitches_top.csv")
liana = read("${WORK_ROOT}/luad_figures/fig5/data/fig5g_liana_lr_pairs.csv")
liana_full = read("${WORK_ROOT}/luad_figures/fig5/data/fig5g_liana_focus_all_senders.csv")
commot = read("${DATA_ROOT}/ST/results/r_data/commot_per_sample_summary.csv")

# Fix cohort label inside COMMOT table if any
if commot is not None:
    for c in ("cohort", "dataset"):
        if c in commot.columns:
            commot[c] = commot[c].replace({"Okamura": "Takano 2024",
                                             "EMTAB13530": "E-MTAB-13530"})
            break

with pd.ExcelWriter(OUT / "Supplementary_Table_S11_TF_LR_COMMOT.xlsx") as xw:
    if tf_top is not None:
        write_titled(xw, "Top state-specific TFs",
            "Table S11a. Top-15 state-specific transcription factors per MP (Fig 3A).",
            "SCENIC AUCell z-scored across MPs; top-15 = highest |z| per MP.",
            tf_top)
    if tf_z is not None:
        write_titled(xw, "TF activity z-score",
            "Table S11b. TF activity z-score matrix (TF x MP1-4) - Fig 3A heatmap data.",
            "Each TF z-scored across 4 MPs; positive = higher in that MP.",
            tf_z)
    if gs is not None:
        write_titled(xw, "GeneSwitches results",
            "Table S11c. GeneSwitches outcome - switch genes along Monocle3 pseudotime (Fig 3C).",
            ("switch_pseudotime_rank = order of switch event along pseudotime; "
             "r2 / pval = logistic-regression fit; direction = up/down."),
            gs)
    if gs_top is not None:
        write_titled(xw, "GeneSwitches top",
            "Table S11d. Top GeneSwitches outputs displayed in Fig 3C.",
            "Filtered to high-confidence binary switches with r2 > 0.3.",
            gs_top)
    if liana is not None:
        write_titled(xw, "LIANA LR pairs",
            "Table S11e. LIANA top ligand-receptor pairs (Fig 5G).",
            ("Aggregated rank from LIANA consensus; lower = stronger inferred "
             "interaction. Pairs limited to Neu_Inflammatory as sender for "
             "clarity; the full sender x LR table is in the next sheet."),
            liana)
    if liana_full is not None:
        write_titled(xw, "LIANA all senders (focus)",
            "Table S11f. LIANA full sender x receiver x LR table for the focus interaction set (Fig 5G support).",
            "All senders considered; subset to focus LRs (EMT-related ligands flagged in Fig 5D).",
            liana_full)
    if commot is not None:
        write_titled(xw, "COMMOT per-section",
            "Table S11g. COMMOT spatial communication summary per section (Fig 7E support).",
            ("Sender / receiver scores aggregated per pathway (OSM, IL1, etc.) "
             "per spatial section. Cohort label uses 'Takano 2024' throughout."),
            commot)

print("\nDONE  -  all tables rebuilt under", OUT)
for f in sorted(OUT.iterdir()):
    if f.name.startswith("Supplementary_Table_S"):
        print(f"  {f.name}  {f.stat().st_size / 1024:.1f} KB")
