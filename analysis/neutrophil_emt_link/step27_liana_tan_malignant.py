"""step27: LIANA cell-cell communication — TAN/NAN → malignant cells stratified by dominant_MP.

Tests two-axis model:
  Axis 1 (EMT-driver):  TAN-1, NAN-2 → Mal_MP1 (immune/EMT infiltrate)
                        focus LR: IL1B-IL1R1/2, CXCL8-CXCR1/2, PLAU-PLAUR,
                                  OSM-OSMR/LIFR, TGFB1-TGFBR1/2, TNF-TNFRSF1A
  Axis 2 (Met-tropic):  TAN-4 → Mal_MP2 (proliferation)
                        focus LR: VEGFA-KDR/NRP1, CCL3/4-CCR1/5

Inputs:
  data/processed/luad_neutrophil_own_raw.h5ad  +  obs from annotated (drop Neutrophils-unclass)
  data/processed/luad_malignant_scored.h5ad

Outputs:
  results/step27_liana_full.csv.gz                  (full LR table)
  results/step27_liana_tan1_to_malignant.csv        (top-30 TAN-1 → Mal_MP*)
  results/step27_liana_tan4_to_malignant.csv        (top-30 TAN-4 → Mal_MP*)
  results/step27_liana_focus_pathways.csv           (only the EMT focus list)
  figures/step27_liana_dotplot_tan1.pdf
  figures/step27_liana_dotplot_tan4.pdf
  figures/step27_liana_focus_heatmap.pdf            (LR × sender-receiver pair)
"""
import os, time, gc
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

t0 = time.time()
DATA = "${PROJECT_ROOT}/data/processed"
RES = "${PROJECT_ROOT}/results"
FIG = f"{RES}/figures"
os.makedirs(FIG, exist_ok=True)

# ----- 1. Build joint anndata -----
print("[1] load neutrophil (annotated obs) + malignant")
neu_ann = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad")
neu_raw = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_raw.h5ad")
# transfer scanvi_predicted to raw
neu_raw.obs["scanvi_predicted"] = neu_ann.obs.loc[neu_raw.obs.index, "scanvi_predicted"].astype(str)
# DROP Neutrophils-unclassified per user instruction
keep = neu_raw.obs["scanvi_predicted"].isin(
    ["TAN-1", "TAN-2", "TAN-3", "TAN-4", "NAN-1", "NAN-2", "NAN-3"])
neu_sub = neu_raw[keep].copy()
print(f"  neu after dropping unclassified: {neu_sub.shape}")
print(f"  per-label counts:\n{neu_sub.obs['scanvi_predicted'].value_counts().to_string()}")

mal = sc.read_h5ad(f"{DATA}/luad_malignant_scored.h5ad")
print(f"  malignant: {mal.shape}")
# build receiver label
mal.obs["dominant_MP"] = mal.obs["dominant_MP"].astype(str)
print(f"  dominant_MP counts:\n{mal.obs['dominant_MP'].value_counts().to_string()}")

# subset malignant to per-MP groups (only keep clear MP labels, drop NaN/other)
mp_keep = mal.obs["dominant_MP"].isin(["MP1", "MP2", "MP3", "MP4"])
mal_sub = mal[mp_keep].copy()
print(f"  malignant after MP filter: {mal_sub.shape}")

# ----- 2. Harmonize var_names + concat -----
print("\n[2] harmonize + concat")
shared = sorted(set(neu_sub.var_names) & set(mal_sub.var_names))
print(f"  shared genes: {len(shared)}")
neu_sub = neu_sub[:, shared].copy()
mal_sub = mal_sub[:, shared].copy()

# unify label
neu_sub.obs["cellgroup"] = "Neu_" + neu_sub.obs["scanvi_predicted"].astype(str)
mal_sub.obs["cellgroup"] = "Mal_" + mal_sub.obs["dominant_MP"].astype(str)
neu_sub.obs["compartment"] = "neutrophil"
mal_sub.obs["compartment"] = "malignant"

# minimal obs to keep before concat
keep_obs = ["cellgroup", "compartment", "sample_id", "dataset", "tissue_type"]
for a in [neu_sub, mal_sub]:
    a.obs = a.obs[[c for c in keep_obs if c in a.obs.columns]].copy()

# raw counts → .X. Both already have raw counts.
import scipy.sparse as sp
def ensure_int(adata):
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    if sp.issparse(adata.X):
        adata.X.data = np.rint(adata.X.data).astype(np.float32)
    else:
        adata.X = np.rint(adata.X).astype(np.float32)
    adata.layers.clear()
ensure_int(neu_sub); ensure_int(mal_sub)

joint = ad.concat([neu_sub, mal_sub], axis=0, join="outer", merge="same",
                  index_unique=None, label="src")
del neu_sub, mal_sub, neu_raw, neu_ann; gc.collect()
print(f"  joint: {joint.shape}")
print(f"  joint cellgroup counts:\n{joint.obs['cellgroup'].value_counts().to_string()}")

# need lognorm in .X for LIANA (it expects normalized data by default; raw counts in layers['counts'])
joint.layers["counts"] = joint.X.copy()
sc.pp.normalize_total(joint, target_sum=1e4)
sc.pp.log1p(joint)

# ----- 3. LIANA -----
print("\n[3] LIANA rank_aggregate")
import liana as li
print(f"  liana version: {li.__version__}")

# rank_aggregate combines 5 methods (CellPhoneDB, NATMI, Connectome, logfc, sca, cellchat-like)
# resource: 'consensus' is the recommended union of curated DBs
li.mt.rank_aggregate(
    joint,
    groupby="cellgroup",
    resource_name="consensus",
    expr_prop=0.1,   # a ligand/receptor must be expressed in >=10% of cells in its group
    verbose=True,
    use_raw=False,
    de_method="wilcoxon",
    n_perms=100,
    seed=0,
    return_all_lrs=False,
)
res = joint.uns["liana_res"]
print(f"  LR result shape: {res.shape}")
print(f"  cols: {list(res.columns)[:20]}")

# save full
res.to_csv(f"{RES}/step27_liana_full.csv.gz", compression="gzip", index=False)

# columns of interest depend on liana version; common: source, target, ligand_complex, receptor_complex,
# magnitude_rank, specificity_rank, lr_means, cellphone_pvals
# normalize column names so we can sort
key_score = "magnitude_rank" if "magnitude_rank" in res.columns else "lr_means"
key_spec = "specificity_rank" if "specificity_rank" in res.columns else "cellphone_pvals"

print(f"  ranking column: {key_score} (lower is stronger by liana convention)")

# ----- 4. focus on TAN-* → Mal_MP* directional pairs -----
print("\n[4] direction filters + focused LR pathways")

def is_neu_to_mal(row):
    return str(row["source"]).startswith("Neu_") and str(row["target"]).startswith("Mal_")

res_n2m = res[res.apply(is_neu_to_mal, axis=1)].copy()
print(f"  Neu→Mal LR rows: {len(res_n2m)}")

# top-30 from TAN-1 / TAN-4 to all Mal_MP*
def top_n(df, src, n=30):
    sub = df[df["source"] == src].sort_values(key_score).head(n)
    return sub

for src in ["Neu_TAN-1", "Neu_TAN-4", "Neu_NAN-2", "Neu_TAN-3", "Neu_NAN-1"]:
    sub = top_n(res_n2m, src, 30)
    safe = src.replace("Neu_", "").replace("-", "")
    sub.to_csv(f"{RES}/step27_liana_top30_{safe}_to_malignant.csv", index=False)
    print(f"  top-30 {src} → Mal: {len(sub)} rows  →  step27_liana_top30_{safe}_to_malignant.csv")

# pathway focus
FOCUS = {
    # Axis 1 (EMT-driver)
    "IL1":     [("IL1B", ["IL1R1", "IL1R2", "IL1RAP"]), ("IL1A", ["IL1R1", "IL1R2"])],
    "CXCL8":   [("CXCL8", ["CXCR1", "CXCR2"])],
    "PLAU":    [("PLAU", ["PLAUR", "LRP1"])],
    "OSM":     [("OSM", ["OSMR", "LIFR", "IL6ST"])],
    "TGFb":    [("TGFB1", ["TGFBR1", "TGFBR2", "TGFBR3"])],
    "TNF":     [("TNF", ["TNFRSF1A", "TNFRSF1B"])],
    # Axis 2 (Met-tropic)
    "VEGFA":   [("VEGFA", ["KDR", "FLT1", "FLT4", "NRP1", "NRP2"])],
    "CCL":     [("CCL3", ["CCR1", "CCR5"]), ("CCL4", ["CCR5"]),
                ("CCL2", ["CCR2"]), ("CCL5", ["CCR5", "CCR1", "CCR3"])],
    "CXCL_other": [("CXCL1", ["CXCR2"]), ("CXCL2", ["CXCR2"])],
    "Other":   [("SPP1", ["CD44", "ITGAV", "ITGB1", "ITGA5"]),
                ("FN1", ["ITGA5", "ITGB1", "CD44"]),
                ("AREG", ["EGFR"]), ("EREG", ["EGFR"])],
}

# helper: detect ligand/receptor cols (LIANA may use *_complex with comma-sep complexes)
lig_col = "ligand_complex" if "ligand_complex" in res.columns else "ligand"
rec_col = "receptor_complex" if "receptor_complex" in res.columns else "receptor"

focus_rows = []
for pw, pairs in FOCUS.items():
    for lig, recs in pairs:
        sub = res_n2m[
            res_n2m[lig_col].astype(str).str.contains(rf"\b{lig}\b", regex=True) &
            res_n2m[rec_col].astype(str).apply(lambda s: any(r in s for r in recs))
        ].copy()
        sub["pathway"] = pw
        focus_rows.append(sub)
focus = pd.concat(focus_rows, ignore_index=True) if focus_rows else pd.DataFrame()
focus = focus.sort_values(["pathway", key_score])
focus.to_csv(f"{RES}/step27_liana_focus_pathways.csv", index=False)
print(f"  focus LR rows: {len(focus)}; pathways: {focus['pathway'].nunique() if len(focus) else 0}")

# ----- 5. visualizations -----
print("\n[5] dotplots + heatmap")
def plot_dot(df, title, outpath, n=25):
    if len(df) == 0:
        print(f"  skip {title}: empty"); return
    df = df.head(n).copy()
    df["lr"] = df[lig_col].astype(str) + " → " + df[rec_col].astype(str)
    df["target_short"] = df["target"].str.replace("Mal_", "")
    pivot_score = df.pivot_table(index="lr", columns="target_short",
                                 values=key_score, aggfunc="min")
    pivot_size = df.pivot_table(index="lr", columns="target_short",
                                values="lr_means" if "lr_means" in df.columns else key_score,
                                aggfunc="mean")
    fig, ax = plt.subplots(figsize=(0.6 * pivot_score.shape[1] + 4,
                                    0.32 * pivot_score.shape[0] + 1.5))
    sns.heatmap(pivot_score, cmap="viridis_r", ax=ax,
                cbar_kws={"label": f"{key_score} (lower = stronger)"})
    ax.set_title(title); ax.set_xlabel("target Mal_MP"); ax.set_ylabel("LR pair")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

# TAN-1 dot
sub1 = res_n2m[res_n2m["source"] == "Neu_TAN-1"].sort_values(key_score)
plot_dot(sub1, "TAN-1 → Malignant top LR pairs", f"{FIG}/step27_liana_dotplot_tan1.pdf")
sub4 = res_n2m[res_n2m["source"] == "Neu_TAN-4"].sort_values(key_score)
plot_dot(sub4, "TAN-4 → Malignant top LR pairs", f"{FIG}/step27_liana_dotplot_tan4.pdf")
subN2 = res_n2m[res_n2m["source"] == "Neu_NAN-2"].sort_values(key_score)
plot_dot(subN2, "NAN-2 → Malignant top LR pairs", f"{FIG}/step27_liana_dotplot_nan2.pdf")

# focus heatmap: rows = LR pair (pathway-grouped); cols = sender→target combo
if len(focus):
    focus["lr"] = focus[lig_col].astype(str) + " → " + focus[rec_col].astype(str)
    focus["sender_target"] = focus["source"].str.replace("Neu_", "") + " → " + focus["target"].str.replace("Mal_", "")
    pivot_focus = focus.pivot_table(index="lr", columns="sender_target",
                                    values=key_score, aggfunc="min")
    fig, ax = plt.subplots(figsize=(0.5 * pivot_focus.shape[1] + 5,
                                    0.30 * pivot_focus.shape[0] + 2))
    sns.heatmap(pivot_focus, cmap="viridis_r", ax=ax,
                cbar_kws={"label": f"{key_score}"})
    ax.set_title("Focus pathways: sender × target LR strength")
    fig.tight_layout()
    fig.savefig(f"{FIG}/step27_liana_focus_heatmap.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/step27_liana_focus_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

# ----- 6. summary print -----
print("\n[6] summary — focus pathway top LR by sender")
if len(focus):
    summary = focus.groupby(["pathway", "source"]).apply(
        lambda d: d.sort_values(key_score).head(1)
    ).reset_index(drop=True)
    print(summary[["pathway", "source", "target", lig_col, rec_col, key_score]].to_string(index=False))

print(f"\nelapsed: {(time.time()-t0)/60:.1f} min")
print("DONE.")
