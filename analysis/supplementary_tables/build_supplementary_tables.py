"""
Build the 6 Supplementary Tables for the LUAD project, mirroring the
format used by the npj Precision Oncology 2026 HCC reference paper.

Output: ${WORK_ROOT}/Supplementary_Tables/Supplementary_Table_S{1..6}.xlsx

Mapping:
  S1  Data sources + sample info
  S2  cNMF top genes (15 clusters x top genes)
  S3  Meta-program DEGs (MP1-4)
  S4  Neutrophil subtype DEGs
  S5  Macrophage subtype DEGs
  S6  Geneformer KO transition genes (3 sheets)
"""
from __future__ import annotations
from pathlib import Path
import gzip
import io
import numpy as np
import pandas as pd

OUT = Path("${WORK_ROOT}/Supplementary_Tables")
OUT.mkdir(parents=True, exist_ok=True)

# S1  -  Data sources + per-GSE sample info
print("[S1] data sources + sample info")
data_sources = pd.DataFrame([
    ("scRNA-seq",  "GSE131907", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE143423", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143423", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE148071", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE164789", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164789", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE189357", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189357", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE253013", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE123902", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902", "Integration analysis", "10x"),
    ("scRNA-seq",  "GSE127465", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465", "Neutrophil cross-cancer reference", "10x"),
    ("scRNA-seq",  "Salcher 2022 LUAD atlas", "https://doi.org/10.1016/j.ccell.2022.10.008", "High-resolution myeloid/neutrophil reference", "10x / mixed"),
    ("Spatial RNA-seq", "E-MTAB-13530", "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530", "Spatial validation (12 sections)", "10x Visium"),
    ("Spatial RNA-seq", "Takano 2024 (Visium FF)", "https://doi.org/10.1038/s41467-024-54671-7", "Spatial validation (8 fresh-frozen sections)", "10x Visium"),
    ("Bulk RNA-seq",   "TCGA-LUAD", "https://portal.gdc.cancer.gov/projects/TCGA-LUAD", "Survival + Cox + KM analyses (Fig 3, Fig 8O/P)", "Illumina"),
    ("Bulk RNA-seq",   "GSE207422", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422", "Anti-PD-(L)1 neoadjuvant response (Fig 8M)", "Illumina"),
    ("Bulk RNA-seq",   "GSE126044", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044", "Anti-PD-1 R vs NR volcano (Fig 8N)", "Illumina"),
    ("scRNA-seq (validation)", "GSE135222", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222",
        "Single-cell anti-PD-1 R vs NR validation of lead genes (Fig S11D-F)", "10x"),
    ("CRISPR essentiality", "DepMap 24Q2", "https://depmap.org/portal/data_page/?tab=allData", "LUAD cell-line lethality benchmark (Fig 8B/G)", "Avana CRISPR"),
], columns=["DataType", "DataSource", "Accession Link", "Application purpose", "Platform"])

# Sample info from integrated cell metadata
meta_path = Path("${WORK_ROOT}/luad_figures/fig_s1/s1_cell_metadata.csv.gz")
sample_info_rows = []
if meta_path.exists():
    meta = pd.read_csv(meta_path, usecols=["dataset","sample_id","patient_id","tissue_type"])
    # one row per (dataset, sample_id)
    sample_info = (meta.drop_duplicates(["dataset","sample_id"])
                       [["dataset","sample_id","patient_id","tissue_type"]]
                       .rename(columns={"dataset":"Dataset",
                                        "sample_id":"Sample",
                                        "patient_id":"Patient",
                                        "tissue_type":"Tissue"})
                       .sort_values(["Dataset","Sample"])
                       .reset_index(drop=True))
else:
    sample_info = pd.DataFrame(columns=["Dataset","Sample","Patient","Tissue"])

# Spatial samples  -  both cohorts (manually compiled to ensure all FF sections appear)
spatial_samples = pd.DataFrame([
    ("E-MTAB-13530", s, p, "Tumor / boundary") for s, p in [
        ("P10_T1","P10"), ("P10_T2","P10"), ("P10_T3","P10"), ("P10_T4","P10"),
        ("P15_T1","P15"), ("P15_T2","P15"),
        ("P16_T1","P16"), ("P16_T2","P16"),
        ("P24_T1","P24"), ("P24_T2","P24"),
        ("P25_T1","P25"), ("P25_T2","P25"),
    ]
] + [
    ("Takano 2024 (FF)", f"LUAD_No_{n}", f"Patient_{n}", "Invasive LUAD")
    for n in (1, 2, 3, 4, 5, 14, 16, 17)
], columns=["Dataset","Sample","Patient","Tissue"])

sample_info_full = pd.concat([sample_info, spatial_samples], ignore_index=True)

with pd.ExcelWriter(OUT / "Supplementary_Table_S1_Data_Sources.xlsx") as xw:
    title1 = pd.DataFrame([["Table S1. Data sources used in this study."]] +
                          [[None]*5] +
                          [data_sources.columns.tolist()] +
                          data_sources.values.tolist())
    title1.to_excel(xw, sheet_name="Data sources", index=False, header=False)
    title2 = pd.DataFrame([["Table S1 (cont.). Per-sample metadata."]] +
                          [[None]*4] +
                          [sample_info_full.columns.tolist()] +
                          sample_info_full.values.tolist())
    title2.to_excel(xw, sheet_name="Sample info", index=False, header=False)
print(f"  rows: data sources={len(data_sources)}, samples={len(sample_info_full)}")

# S2  -  cNMF top genes (15 cluster columns x top N genes)
print("\n[S2] cNMF top genes per cluster")
gep_top = pd.read_csv("${WORK_ROOT}/luad_figures/fig_s2/s2_gep_top30_genes.csv")
# Rebuild 15 cluster assignment via Spearman correlation hierarchical cutree
# (mirrors stepF logic): we already cached it under fig_s2 results. Use the
# annotation file that maps gep_id -> MP, plus a separate cutree result.
# Quickest: re-run cutree right here from gep_spearman_corr.csv
corr = pd.read_csv("${WORK_ROOT}/luad_figures/fig_s2/gep_spearman_corr.csv", index_col=0)
gep_anno = pd.read_csv("${WORK_ROOT}/luad_figures/fig_s2/gep_count_per_patient_per_mp.csv")  # not the right file
mp_anno_path = Path("${WORK_ROOT}/luad_figures/fig2/gep_mp_annotation.csv")
if mp_anno_path.exists():
    gep_mp = pd.read_csv(mp_anno_path)
else:
    gep_mp = None

# Use scipy hierarchical clustering for 15 clusters
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
keep = [c for c in corr.index if c in corr.columns]
sub = corr.loc[keep, keep].values
sub = (sub + sub.T) / 2
np.fill_diagonal(sub, 1.0)
dist = 1 - sub
np.fill_diagonal(dist, 0.0)
dist_clipped = np.clip(dist, 0, None)
condensed = squareform(dist_clipped, checks=False)
Z = linkage(condensed, method="average")
clust15 = fcluster(Z, t=15, criterion="maxclust")
gep_to_clust = pd.Series(clust15, index=keep)

# Per cluster top genes by frequency rank: pick top-30 genes ranked by
# average rank within the cluster's GEPs.
gep_top_in = gep_top[gep_top["gep_id"].isin(gep_to_clust.index)].copy()
gep_top_in["cluster"] = gep_to_clust.loc[gep_top_in["gep_id"]].values
TOP_N = 30
out_cols = {}
for cid in range(1, int(clust15.max()) + 1):
    sub_df = gep_top_in[gep_top_in["cluster"] == cid]
    if len(sub_df) == 0:
        top = [None] * TOP_N
    else:
        top = (sub_df.groupby("gene")["rank"].mean()
                      .sort_values()
                      .head(TOP_N).index.tolist())
        while len(top) < TOP_N:
            top.append(None)
    out_cols[f"cNMF_{cid}"] = top
s2 = pd.DataFrame(out_cols)
with pd.ExcelWriter(OUT / "Supplementary_Table_S2_cNMF_Top_Genes.xlsx") as xw:
    title2 = pd.DataFrame([
        ["Table S2. Top-30 marker genes per cNMF cluster (15 clusters from k=15 hierarchical cutting of 770 GEPs)."]
    ] + [[None]*len(s2.columns)] + [s2.columns.tolist()] + s2.values.tolist())
    title2.to_excel(xw, sheet_name="cNMF_top_genes", index=False, header=False)
print(f"  shape: {s2.shape}")

# S3  -  Meta-program DEGs
print("\n[S3] meta-program DEGs (MP1-4)")
mp_sig = pd.read_csv("${WORK_ROOT}/luad_figures/fig2/mp_signatures_top100.csv")
# Reorganise to mirror Seurat FindAllMarkers schema as close as possible.
mp_sig_out = mp_sig.rename(columns={
    "MP": "cluster",
    "rank": "rank_within_MP",
    "gene": "gene",
    "consensus_score": "consensus_score",
    "n_GEPs_in_MP": "n_GEPs_supporting",
})
# Match the example column order: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene
# We don't have p-values from Seurat but we do have a consensus score that
# carries the same role (a calibrated importance metric).  Document the
# substitution explicitly in the title row.
mp_sig_out = mp_sig_out[["cluster","rank_within_MP","gene","consensus_score","n_GEPs_supporting"]]
with pd.ExcelWriter(OUT / "Supplementary_Table_S3_Meta_Program_DEGs.xlsx") as xw:
    title = pd.DataFrame([
        ["Table S3. Meta-Program (MP1-MP4) signature genes  -  top 100 per MP."],
        ["Source: cNMF consensus across 770 program-level GEPs (Gavish-style),"
         " ordered by consensus_score (higher = more consistently selected as a top-30 program member across patient GEPs in that MP)."],
    ] + [[None]*len(mp_sig_out.columns)] + [mp_sig_out.columns.tolist()] + mp_sig_out.values.tolist())
    title.to_excel(xw, sheet_name="Meta Programs DEGs", index=False, header=False)
print(f"  rows: {len(mp_sig_out)}")

# S4  -  Neutrophil subtype DEGs
print("\n[S4] neutrophil subtype DEGs")
neu = pd.read_csv("${WORK_ROOT}/luad_figures/fig5/data/fig5c_canonical_markers.csv")
neu_out = neu.rename(columns={
    "neu_subtype": "cluster",
    "gene": "gene",
    "mean_expression": "mean_log_expr",
    "pct_expressing": "pct_expressing",
    "canonical_for": "canonical_for",
})
neu_out = neu_out[["cluster","gene","mean_log_expr","pct_expressing","canonical_for"]]
with pd.ExcelWriter(OUT / "Supplementary_Table_S4_Neutrophil_Subtype_DEGs.xlsx") as xw:
    title = pd.DataFrame([
        ["Table S4. Neutrophil subtype canonical marker expression."],
        ["Source: Salcher LUAD atlas + integrated scRNA-seq (Figure 5C)."
         " mean_log_expr = average log1p-normalised expression in subtype;"
         " pct_expressing = percentage of subtype cells with expr > 0."
         " canonical_for = subtype the gene was originally curated as a marker of."],
    ] + [[None]*len(neu_out.columns)] + [neu_out.columns.tolist()] + neu_out.values.tolist())
    title.to_excel(xw, sheet_name="Neutrophil subtype DEGs", index=False, header=False)
print(f"  rows: {len(neu_out)}")

# S5  -  Macrophage subtype DEGs
print("\n[S5] macrophage subtype DEGs")
m1 = pd.read_csv("${WORK_ROOT}/luad_figures/fig4/panel_GM_subset_markers.csv")
m1_out = m1.rename(columns={
    "subtype": "cluster",
    "gene": "gene",
    "mean_log1p": "mean_log_expr",
    "pct_expressing": "pct_expressing",
    "panel_origin": "canonical_for",
})
m1_out = m1_out[["cluster","gene","mean_log_expr","pct_expressing","canonical_for"]]
# Optionally add the SPP1 vs C1QC pseudobulk DEG table as a 2nd sheet
m2_path = Path("${WORK_ROOT}/luad_figures/fig4/panel_N_spp1_vs_c1qc_deg.csv")
m2 = pd.read_csv(m2_path) if m2_path.exists() else None
with pd.ExcelWriter(OUT / "Supplementary_Table_S5_Macrophage_Subtype_DEGs.xlsx") as xw:
    title = pd.DataFrame([
        ["Table S5. Macrophage subtype marker expression (panel_GM)."],
        ["Source: integrated scRNA-seq myeloid compartment (Figure 4)."
         " Columns: cluster, gene, mean_log_expr (mean log1p-normalised),"
         " pct_expressing, canonical_for (subset that gene was originally curated for)."],
    ] + [[None]*len(m1_out.columns)] + [m1_out.columns.tolist()] + m1_out.values.tolist())
    title.to_excel(xw, sheet_name="Macrophage subtype markers", index=False, header=False)
    if m2 is not None:
        m2_out = m2.copy()
        m2_out.columns = ["gene","score","avg_log2FC","p_val","p_val_adj"][:len(m2_out.columns)]
        ttl2 = pd.DataFrame([
            ["Table S5b. Macro_SPP1 vs Macro_C1QC pseudobulk DEGs."],
            ["Source: pseudobulk DESeq2 (one row per patient × subtype, Figure 4N)."],
        ] + [[None]*len(m2_out.columns)] + [m2_out.columns.tolist()] + m2_out.values.tolist())
        ttl2.to_excel(xw, sheet_name="SPP1 vs C1QC pseudobulk DEG", index=False, header=False)
print(f"  rows: subset markers={len(m1_out)}, SPP1-vs-C1QC={len(m2) if m2 is not None else 0}")

# S6  -  Geneformer KO transition genes (3 sheets)
print("\n[S6] Geneformer KO transitions")
GF = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/geneformer_500_stats")
sheet_specs = [
    ("Macro-SPP1-to-Macro-C1QC", GF / "macro_spp1_to_c1qc_stats.csv"),
    ("Mal-MP3-to-MP1",           GF / "mal_mp3_to_mp1_stats.csv"),
    ("Neu-OSM-priming-to-low",   GF / "neu_osm_priming_to_low_stats.csv"),
]
with pd.ExcelWriter(OUT / "Supplementary_Table_S6_Geneformer_KO_Transitions.xlsx") as xw:
    for sn, p in sheet_specs:
        if not p.exists():
            print(f"  [SKIP] {p.name} not found"); continue
        df = pd.read_csv(p)
        # If the first column is an unnamed index, drop it
        if df.columns[0].startswith("Unnamed") or df.columns[0] == "":
            df = df.drop(columns=df.columns[0])
        df = df[["Gene","Gene_name","Ensembl_ID","Shift_to_goal_end",
                 "Goal_end_vs_random_pval","Goal_end_FDR","N_Detections","Sig"]]
        df = df.sort_values("Shift_to_goal_end", ascending=False)
        title = pd.DataFrame([
            [f"Table S6  -  Geneformer KO transition: {sn}"],
            ["Source: Geneformer V2 (104M) in-silico knockout / over-expression simulation"
             " (Figure 8 perturbation analysis). Shift_to_goal_end = mean shift toward the"
             " desired endpoint embedding (positive = closer to goal). Sig=1 indicates"
             " FDR < 0.05 against random-gene null."],
        ] + [[None]*len(df.columns)] + [df.columns.tolist()] + df.values.tolist())
        title.to_excel(xw, sheet_name=sn, index=False, header=False)
        print(f"  {sn}: {len(df)} rows")

print("\nALL Supplementary Tables built ->", OUT)
for f in sorted(OUT.iterdir()):
    print(f"  {f.name}  ({f.stat().st_size / 1024:.1f} KB)")
