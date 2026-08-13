"""Build Supplementary Tables S7-S11 (Tier 1 + Tier 2) for the LUAD project."""
from __future__ import annotations
from pathlib import Path
import pandas as pd

OUT = Path("${WORK_ROOT}/Supplementary_Tables")
OUT.mkdir(parents=True, exist_ok=True)


def write_titled_sheet(xw, sheet, title, subtitle, df):
    """Write df with a 1-row title, a 1-row subtitle, blank, header, data."""
    rows = [[title] + [None]*(len(df.columns)-1)]
    if subtitle:
        rows.append([subtitle] + [None]*(len(df.columns)-1))
    rows.append([None]*len(df.columns))
    rows.append(list(df.columns))
    rows.extend(df.values.tolist())
    pd.DataFrame(rows).to_excel(xw, sheet_name=sheet[:31], index=False, header=False)


def safe_read(p, **kw):
    p = Path(p)
    if not p.exists():
        print(f"  [SKIP] missing {p}")
        return None
    return pd.read_csv(p, **kw)


# S7 Survival statistics
print("[S7] Survival statistics")
fig3 = Path("${WORK_ROOT}/luad_figures/fig3")
fig8d = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data")

cox_mv = safe_read(fig3 / "tcga_luad_mp_cox_multivariate.csv")
cox_uv = safe_read(fig3 / "tcga_luad_mp_cox_univariate.csv")
km_med = safe_read(fig3 / "tcga_luad_mp_km_logrank.csv")
km_opt = safe_read(fig3 / "tcga_optimal_cutoff_km.csv")
km_risk = safe_read(fig3 / "tcga_risk_score_km.csv")
km_srsf9 = safe_read(fig8d / "8O_km_SRSF9_stats.csv")
km_sec61g = safe_read(fig8d / "8P_km_SEC61G_stats.csv")

with pd.ExcelWriter(OUT / "Supplementary_Table_S7_Survival_Statistics.xlsx") as xw:
    if cox_uv is not None:
        write_titled_sheet(xw, "Cox univariate",
            "Table S7a. TCGA-LUAD Cox proportional-hazards univariate regression for each MP score (Fig 3C input).",
            "Source: TCGA-LUAD overall survival, MP1-4 scored on bulk RNA-seq via single-sample GSEA.",
            cox_uv)
    if cox_mv is not None:
        write_titled_sheet(xw, "Cox multivariate",
            "Table S7b. TCGA-LUAD Cox multivariate regression (Age + Stage + MP1-4)   Fig 3C forest plot.",
            "Source: TCGA-LUAD; covariates Age (continuous), Stage (numeric I-IV), MP1-MP4 z-scored.",
            cox_mv)
    if km_med is not None:
        write_titled_sheet(xw, "KM median split",
            "Table S7c. TCGA-LUAD KM log-rank by median MP score split (Fig 3B).",
            "n_high / n_low / n_events_high / n_events_low / log-rank p / HR (high vs low).",
            km_med)
    if km_opt is not None:
        write_titled_sheet(xw, "KM optimal cutoff",
            "Table S7d. TCGA-LUAD KM with optimal-cutoff (max-rank) split per MP (Fig 3D).",
            "Optimal cutoff selected to maximise log-rank statistic (constrained to 30-70 percentile).",
            km_opt)
    if km_risk is not None:
        write_titled_sheet(xw, "Risk score KM",
            "Table S7e. TCGA-LUAD risk-score (linear combination of MP coefficients) KM (Fig 3F).",
            "Risk score = beta * MP_score combination from multivariate Cox; tertile split.",
            km_risk)
    # Lead gene KM (Fig 8O/P)
    pieces = []
    for stat, gene in [(km_sec61g, "SEC61G"), (km_srsf9, "SRSF9")]:
        if stat is None: continue
        s = stat.copy()
        if "gene" in s.columns:
            s = s.drop(columns=["gene"])
        s.insert(0, "gene", gene)
        pieces.append(s)
    if pieces:
        lg = pd.concat(pieces, ignore_index=True)
        write_titled_sheet(xw, "Lead gene KM",
            "Table S7f. TCGA-LUAD KM for the two lead genes SEC61G and SRSF9 (Fig 8O/P).",
            "High vs Low split via optimal cutoff; covariate-adjusted log-rank p reported.",
            lg)
print("  OK")


# S8 Spatial validation statistics
print("[S8] Spatial validation statistics")
r_data = Path("${DATA_ROOT}/ST/results/r_data")
roi_agg = safe_read(r_data / "roi_vs_nonroi_aggregate_pvalues.csv")
roi_sec = safe_read(r_data / "roi_vs_nonroi_stats_with_pvalues.csv")
roi_consist = safe_read(r_data / "roi_per_section_consistency.csv")
misty = safe_read(r_data / "misty_aggregated_importance.csv")
sec_mp = safe_read(fig8d / "8I_section_MP_scores.csv")

# also include the lead-gene per-section table from Fig 8 work
lead_sec = safe_read(fig8d / "per_section_lead_gene_sig.csv")
lead_coh = safe_read(fig8d / "per_cohort_lead_gene_sig.csv")

with pd.ExcelWriter(OUT / "Supplementary_Table_S8_Spatial_Validation_Statistics.xlsx") as xw:
    if roi_agg is not None:
        # Replace internal "Okamura" with "Takano 2024" for display
        if "cohort" in roi_agg.columns:
            roi_agg["cohort"] = roi_agg["cohort"].replace({"Okamura": "Takano 2024",
                                                             "EMTAB13530": "E-MTAB-13530"})
        write_titled_sheet(xw, "ROI vs non-ROI cohort",
            "Table S8a. ROI vs non-ROI per-cohort statistics (Fig 7Q, Fig S10J).",
            "Per-spot Mann-Whitney U test (two-sided) pooled within cohort across sections; BH-FDR within cohort across all 67 metrics.",
            roi_agg)
    if roi_sec is not None:
        if "cohort" in roi_sec.columns:
            roi_sec["cohort"] = roi_sec["cohort"].replace({"Okamura": "Takano 2024",
                                                             "EMTAB13530": "E-MTAB-13530"})
        write_titled_sheet(xw, "ROI vs non-ROI per-section",
            "Table S8b. ROI vs non-ROI per-section statistics (Fig S10K input).",
            "Mann-Whitney U test per (sample, metric); BH-FDR within each section across metrics. Sections with <5 ROI or <5 non-ROI spots dropped.",
            roi_sec)
    if roi_consist is not None:
        write_titled_sheet(xw, "Per-section consistency",
            "Table S8c. Per-metric consistency summary across all sections (Fig S10K caption support).",
            "n_sig_strong = sections with FDR<0.01; n_sig_any = sections with FDR<0.05; sign_consistent = sections agreeing with metric's median delta sign.",
            roi_consist)
    if misty is not None:
        write_titled_sheet(xw, "MISTy importance",
            "Table S8d. MISTy cell-type interaction importances (Fig 7H).",
            "Aggregated importance per (target, predictor) across views (intra/juxta/para); higher = stronger spatial dependency.",
            misty)
    if sec_mp is not None:
        write_titled_sheet(xw, "Per-section MP scores",
            "Table S8e. Per-section mean MP1-4 / dominant MP scores (Fig 8I auxiliary).",
            "Mean spot-level MP score per (cohort, sample); used to assign R-surrogate vs NR-surrogate sections.",
            sec_mp)
    if lead_coh is not None:
        write_titled_sheet(xw, "Lead gene cohort sig",
            "Table S8f. Tumor-intrinsic ROI vs non-ROI per-cohort Mann-Whitney + FDR for the 3 lead genes (Fig S11Q).",
            "ROI = z(Malignant)>0.5 AND z(MP3)>0.5; per-spot test pooled within cohort; BH-FDR within cohort across 3 genes.",
            lead_coh)
    if lead_sec is not None:
        write_titled_sheet(xw, "Lead gene per-section sig",
            "Table S8g. Tumor-intrinsic ROI vs non-ROI per-section Mann-Whitney + FDR for the 3 lead genes (Fig S11R).",
            "Per-section spot-level test; BH-FDR within each section across 3 genes.",
            lead_sec)
print("  OK")


# S9 Treatment / clinical validation
print("[S9] Treatment / clinical validation")
depmap = safe_read(fig8d / "8B_depmap_stats.csv")
expr_eff = safe_read(fig8d / "8G_expr_vs_effect.csv")
tcga_tvn = safe_read(fig8d / "8D_tcga_TvN_stats.csv")
mpr = safe_read(fig8d / "8M_GSE207422_stats.csv")
volcano = safe_read(fig8d / "8N_GSE126044_volcano.csv")
gse135 = safe_read(fig8d / "S11D_GSE135222_stats.csv")

with pd.ExcelWriter(OUT / "Supplementary_Table_S9_Treatment_Clinical_Validation.xlsx") as xw:
    if depmap is not None:
        write_titled_sheet(xw, "DepMap LUAD essentiality",
            "Table S9a. DepMap 24Q2 CRISPR essentiality of lead genes in LUAD cell lines (Fig 8B).",
            "Gene effect (Chronos score, more negative = more essential); LUAD cohort vs non-LUAD; one-sided Mann-Whitney.",
            depmap)
    if expr_eff is not None:
        write_titled_sheet(xw, "Expr × essentiality",
            "Table S9b. Per-cell-line gene expression × CRISPR effect for lead genes (Fig 8G scatter).",
            "log2(TPM+1) vs gene_effect; LUAD lines only; Spearman rho per gene.",
            expr_eff)
    if tcga_tvn is not None:
        write_titled_sheet(xw, "TCGA Tumor vs Normal",
            "Table S9c. TCGA-LUAD Tumor vs Normal differential expression for lead genes (Fig 8D-F).",
            "Wilcoxon two-sided per gene; log2FC = log2(mean Tumor TPM+1) − log2(mean Normal TPM+1).",
            tcga_tvn)
    if mpr is not None:
        write_titled_sheet(xw, "GSE207422 MPR vs NMPR",
            "Table S9d. GSE207422 (neoadjuvant chemo-IO) MPR vs NMPR per-gene statistics (Fig 8M).",
            "Per-gene Wilcoxon two-sided + AUC of MPR vs NMPR responder classification.",
            mpr)
    if volcano is not None:
        write_titled_sheet(xw, "GSE126044 R vs NR",
            "Table S9e. GSE126044 (anti-PD-1) Responder vs Non-Responder full DE (Fig 8N volcano).",
            "Wilcoxon per-gene; log2FC = R − NR; -log10(p) on y-axis. Lead genes flagged via is_lead column.",
            volcano)
    if gse135 is not None:
        write_titled_sheet(xw, "GSE135222 single-cell",
            "Table S9f. GSE135222 single-cell anti-PD-1 R vs NR validation of lead genes (Fig S11D-F).",
            "Per-gene Wilcoxon on log-normalised expression in tumor compartment; R vs NR.",
            gse135)
print("  OK")


# S10 Pathway / hallmark enrichment
print("[S10] Pathway / hallmark enrichment")
fig2 = Path("${WORK_ROOT}/luad_figures/fig2")
hall_nes = safe_read(fig2 / "hallmark_nes_heatmap.csv")
hall_fdr = safe_read(fig2 / "hallmark_fdr_heatmap.csv")
gavish = safe_read(fig2 / "gavish_top_matches.csv")
gavish_overlap = safe_read(fig2 / "gavish_overlap.csv")

# PROGENy per-section means: use both step08 and step09
prog_e = safe_read("${DATA_ROOT}/ST/results/step05_progeny/per_sample_mean_progeny.csv")
prog_o = safe_read("${DATA_ROOT}/ST/results/step09_okamura_validation/per_sample_mean_progeny.csv")
prog_pieces = []
if prog_e is not None:
    prog_e2 = prog_e.copy(); prog_e2.insert(0, "cohort", "E-MTAB-13530"); prog_pieces.append(prog_e2)
if prog_o is not None:
    prog_o2 = prog_o.copy(); prog_o2.insert(0, "cohort", "Takano 2024"); prog_pieces.append(prog_o2)
prog_all = pd.concat(prog_pieces, ignore_index=True) if prog_pieces else None

with pd.ExcelWriter(OUT / "Supplementary_Table_S10_Pathway_Enrichment.xlsx") as xw:
    if hall_nes is not None:
        write_titled_sheet(xw, "Hallmark NES",
            "Table S10a. Hallmark gene-set NES per Meta-Program (Fig 2F heatmap rows).",
            "GSEA NES on MP scores vs all-other; positive = enriched in MP. 50 Hallmark sets × MP1-4.",
            hall_nes)
    if hall_fdr is not None:
        write_titled_sheet(xw, "Hallmark FDR",
            "Table S10b. BH-adjusted p-values matching the NES table (Fig 2F stars overlay).",
            "FDR < 0.05 / 0.01 / 0.001 = * / ** / *** in the figure.",
            hall_fdr)
    if gavish is not None:
        write_titled_sheet(xw, "Gavish pan-cancer alignment",
            "Table S10c. Top Gavish 41-MP pan-cancer match per cNMF cluster (Fig 2I).",
            "Cosine similarity between cNMF cluster mean profile and each Gavish 2023 pan-cancer MP signature.",
            gavish)
    if gavish_overlap is not None:
        write_titled_sheet(xw, "Gavish overlap",
            "Table S10d. Gene-set overlap matrix between cNMF clusters and Gavish 2023 MPs (Fig 2I support).",
            "Overlap = |intersect| / min(|set_A|, |set_B|) using top-30 genes of each program.",
            gavish_overlap)
    if prog_all is not None:
        write_titled_sheet(xw, "PROGENy per-section",
            "Table S10e. Mean PROGENy 14-pathway activity per spatial section (Fig 7I, Fig S10H).",
            "Mean MLM-derived score per section across 14 pathways; both ST cohorts pooled.",
            prog_all)
print("  OK")


# S11 TF / GeneSwitches / Ligand-Receptor / COMMOT
print("[S11] TF / GeneSwitches / LR / COMMOT")
tf_top = safe_read(fig3 / "tf_state_specific_top15.csv")
tf_z = safe_read(fig3 / "tf_activity_mp_zscore.csv")
gs = safe_read(fig3 / "geneswitches_results.csv")
gs_top = safe_read(fig3 / "geneswitches_top.csv")
liana = safe_read("${WORK_ROOT}/luad_figures/fig5/data/fig5g_liana_lr_pairs.csv")
liana_full = safe_read("${WORK_ROOT}/luad_figures/fig5/data/fig5g_liana_focus_all_senders.csv")
commot = safe_read("${DATA_ROOT}/ST/results/r_data/commot_per_sample_summary.csv")

with pd.ExcelWriter(OUT / "Supplementary_Table_S11_TF_LR_COMMOT.xlsx") as xw:
    if tf_top is not None:
        write_titled_sheet(xw, "Top state-specific TFs",
            "Table S11a. Top-15 state-specific transcription factors per MP (Fig 3A).",
            "SCENIC AUCell z-scored across MPs; top-15 = highest |z| per MP.",
            tf_top)
    if tf_z is not None:
        write_titled_sheet(xw, "TF activity z-score",
            "Table S11b. TF activity z-score matrix (TF × MP1-4)   Fig 3A heatmap data.",
            "Each TF z-scored across 4 MPs; positive = higher in that MP.",
            tf_z)
    if gs is not None:
        write_titled_sheet(xw, "GeneSwitches results",
            "Table S11c. GeneSwitches outcome   switch genes along Monocle3 pseudotime (Fig 3 axis).",
            "switch_pseudotime_rank = order of switch event along pseudotime; r2 / pval = logistic-regression fit; direction = up/down.",
            gs)
    if gs_top is not None:
        write_titled_sheet(xw, "GeneSwitches top",
            "Table S11d. Top GeneSwitches outputs displayed in Fig 3C.",
            "Filtered to high-confidence binary switches with r2 > 0.3.",
            gs_top)
    if liana is not None:
        write_titled_sheet(xw, "LIANA LR pairs",
            "Table S11e. LIANA top ligand-receptor pairs (Fig 5G).",
            "Aggregated rank from LIANA consensus; lower = stronger inferred interaction.",
            liana)
    if liana_full is not None:
        write_titled_sheet(xw, "LIANA all senders (focus)",
            "Table S11f. LIANA full sender × receiver × LR table for the focus interaction set (Fig 5G support).",
            "All senders considered; subset to focus LRs (EMT-related ligands flagged in Fig 5D).",
            liana_full)
    if commot is not None:
        write_titled_sheet(xw, "COMMOT per-section",
            "Table S11g. COMMOT spatial communication summary per section (Fig 7E support).",
            "Sender / receiver scores aggregated per pathway (OSM, IL1, etc.) per spatial section.",
            commot)
print("  OK")

print("\nALL Tier 1+2 tables written ->", OUT)
for f in sorted(OUT.iterdir()):
    if f.name.startswith("Supplementary_Table_S"):
        print(f"  {f.name}  ({f.stat().st_size / 1024:.1f} KB)")
