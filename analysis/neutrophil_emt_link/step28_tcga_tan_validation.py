"""step28: TCGA-LUAD external validation of two-axis TAN→tumor model.

Three analyses:
  1. ssGSEA TAN-1/2/3/4, NAN-1/2/3 sigs on TCGA-LUAD TPM, then
     Spearman correlate vs MP1-5 (existing) + Hallmark EMT
  2. Stratified KM:
     - TAN-1_high × MP1_high vs others (4-group)
     - TAN-4_high × MP2_high vs others (4-group)
  3. Multivariate Cox: TAN-1 + MP1 + age + stage + gender

Inputs:
  ${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv
  ${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv
  ${WORK_ROOT}/luad_figures/fig3/tcga_luad_mp_ssgsea.csv.gz   (existing MP scores per sample)
  results/step26_tan_signatures.csv   (TAN/NAN signatures, padj<0.05+lfc>0.5)
  results/step25e_markers_scanvi_fullgene.csv  (full DE for TAN-2 relax)

Outputs:
  results/step28_tcga_tan_ssgsea.csv.gz
  results/step28_tcga_correlation_matrix.csv
  results/step28_km_*_logrank.csv
  results/step28_cox_multivariate.csv
  figures/step28_tcga_correlation_heatmap.pdf
  figures/step28_km_tan1_mp1_4group.pdf
  figures/step28_km_tan4_mp2_4group.pdf
  figures/step28_cox_forest.pdf
"""
import os, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

t0 = time.time()
TPM_CSV = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv")
CLIN_CSV = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
MP_SCORES = Path("${WORK_ROOT}/luad_figures/fig3/tcga_luad_mp_ssgsea.csv.gz")
SIG_CSV = Path("${PROJECT_ROOT}/results/step26_tan_signatures.csv")
SIG_FULL = Path("${PROJECT_ROOT}/results/step25e_markers_scanvi_fullgene.csv")
RES = Path("${PROJECT_ROOT}/results")
FIG = RES / "figures"
FIG.mkdir(exist_ok=True, parents=True)

# ----- 1. Build TAN/NAN signatures (relax TAN-2 to top-50 by score) -----
print("[1] build signatures (TAN-2 relaxed) + apply functional rename")
sig = pd.read_csv(SIG_CSV)
sig_full = pd.read_csv(SIG_FULL)

# Functional rename per user 2026-04-27
RENAME = {
    "TAN-1": "Neu_Inflammatory",
    "TAN-2": "Neu_IFN_response",
    "TAN-3": "Neu_Angiogenic",
    "TAN-4": "Neu_Metastatic",
    "NAN-1": "Neu_ECM_remodeling",
    "NAN-2": "Neu_OSM_priming",
    "NAN-3": "Neu_OSM_low",
}
GROUPS = ["TAN-1", "TAN-2", "TAN-3", "TAN-4", "NAN-1", "NAN-2", "NAN-3"]
ORDER = [RENAME[g] for g in GROUPS]
gene_sets = {}
for g in GROUPS:
    new_name = RENAME[g]
    s = sig[sig["group"] == g]
    if len(s) < 5:
        s2 = sig_full[sig_full["group"] == g].sort_values("score", ascending=False).head(50)
        gene_sets[new_name] = s2["gene"].tolist()
        print(f"  {new_name} (was {g}): relaxed to top-50 by score (strict had {len(s)})")
    else:
        gene_sets[new_name] = s["gene"].tolist()[:50]
        print(f"  {new_name} (was {g}): {len(gene_sets[new_name])} genes")

# Hallmark EMT — pull from gseapy/MSigDB or hardcode top genes
HALLMARK_EMT = """ABI3BP ACTA2 ADAM12 ANPEP APLP1 AREG BASP1 BDNF BGN BMP1 CADM1 CALD1 CALU CAP2 CAPG CAPN1 CAV1 CD44 CD59 CDH11 CDH2 CDH6 COL11A1 COL12A1 COL16A1 COL1A1 COL1A2 COL3A1 COL4A1 COL4A2 COL5A1 COL5A2 COL5A3 COL6A2 COL6A3 COL7A1 COL8A2 COMP COPA CRLF1 CTGF CTHRC1 CXCL1 CXCL12 CXCL6 CXCL8 CYR61 DAB2 DCN DKK1 DPYSL3 ECM1 ECM2 EDIL3 EFEMP2 ELN EMP3 ENO2 FAP FAS FBLN1 FBLN2 FBLN5 FBN1 FBN2 FERMT2 FGF2 FLNA FMOD FN1 FOXC2 FSTL1 FSTL3 FUCA1 FZD8 GADD45A GADD45B GAS1 GEM GJA1 GLIPR1 GPC1 GPX7 GREM1 HTRA1 ID2 IGFBP2 IGFBP3 IGFBP4 IL15 IL32 IL6 INHBA ITGA2 ITGA5 ITGAV ITGB1 ITGB3 ITGB5 JUN LAMA1 LAMA2 LAMA3 LAMC1 LAMC2 LGALS1 LOX LOXL1 LOXL2 LRP1 LRRC15 LUM MAGEE1 MATN2 MATN3 MCM7 MEST MFAP5 MGP MMP1 MMP14 MMP2 MMP3 MSX1 MYL9 MYLK NID2 NNMT NOTCH2 NT5E NTM OXTR PCOLCE PCOLCE2 PDGFRB PDLIM4 PFN2 PLAUR PLOD1 PLOD2 PLOD3 PMEPA1 PMP22 POSTN PPIB PRRX1 PRSS2 PTHLH PTX3 PVR QSOX1 RGS4 RHOB SAT1 SCG2 SDC1 SDC4 SERPINE1 SERPINE2 SERPINH1 SFRP1 SFRP4 SGCB SGCD SGCG SLC6A8 SLIT2 SLIT3 SNAI2 SNTB1 SPARC SPOCK1 SPP1 TAGLN TFPI2 TGFB1 TGFBI TGFBR3 TGM2 THBS1 THBS2 THY1 TIMP1 TIMP3 TNC TNFAIP3 TNFRSF11B TNFRSF12A TPM1 TPM2 TPM4 VCAM1 VCAN VEGFA VEGFC VIM WIPF1 WNT5A""".split()
gene_sets["EMT_Hallmark"] = HALLMARK_EMT

# ----- 2. ssGSEA on TCGA TPM -----
print("\n[2] load TCGA TPM + clinical")
tpm = pd.read_csv(TPM_CSV, index_col=0)
clin = pd.read_csv(CLIN_CSV)
clin_pt = clin[clin["sample_type"] == "Primary Tumor"].copy()
pt_samples = [s for s in tpm.columns if s in set(clin_pt["sample_barcode"])]
tpm = tpm[pt_samples]
if not tpm.index.is_unique:
    tpm = tpm.groupby(tpm.index).max()
print(f"  TPM after dedup: {tpm.shape}")

# gene availability check
print("\n  signature gene availability in TCGA TPM:")
for k, gs in gene_sets.items():
    have = [g for g in gs if g in tpm.index]
    print(f"    {k}: {len(have)}/{len(gs)}")

print("\n[2b] ssGSEA via gseapy")
import gseapy as gp
expr = np.log2(tpm + 1.0).astype("float32")
ss = gp.ssgsea(
    data=expr,
    gene_sets=gene_sets,
    outdir=None,
    sample_norm_method="rank",
    no_plot=True,
    min_size=5,
    max_size=5000,
    permutation_num=0,
    seed=0,
    threads=8,
)
scores = ss.res2d.pivot_table(index="Name", columns="Term", values="NES").astype(float)
scores.index.name = "sample_barcode"
scores.to_csv(RES / "step28_tcga_tan_ssgsea.csv.gz", compression="gzip")
print(f"  TAN/NAN scores shape: {scores.shape}")

# ----- 3. Merge with existing MP scores + correlation matrix -----
print("\n[3] correlate TAN/NAN sigs vs MP1-5 + EMT")
mp_scores = pd.read_csv(MP_SCORES, index_col=0)
print(f"  MP scores shape: {mp_scores.shape}")
combined = scores.join(mp_scores, how="inner")
print(f"  combined shape: {combined.shape}")
combined.to_csv(RES / "step28_tcga_combined_scores.csv.gz", compression="gzip")

# correlation
sig_cols = list(gene_sets.keys())
mp_cols = [c for c in combined.columns if c.startswith("MP")]
target_cols = mp_cols  # don't include EMT_Hallmark vs itself

cor_rows = []
for sc_col in sig_cols:
    for mc in mp_cols:
        rho, p = spearmanr(combined[sc_col], combined[mc])
        cor_rows.append({"sig": sc_col, "MP": mc, "spearman_rho": rho, "p": p, "n": len(combined)})
cor_df = pd.DataFrame(cor_rows)
cor_df.to_csv(RES / "step28_tcga_correlation_matrix.csv", index=False)
rho_mat = cor_df.pivot_table(index="sig", columns="MP", values="spearman_rho")
p_mat = cor_df.pivot_table(index="sig", columns="MP", values="p")
print("\n  ρ matrix:")
print(rho_mat.round(3).to_string())

# heatmap
def annot_star(rho, p):
    s = f"{rho:.2f}"
    if p < 0.001: s += "***"
    elif p < 0.01: s += "**"
    elif p < 0.05: s += "*"
    return s
annot = pd.DataFrame(np.vectorize(annot_star)(rho_mat.values, p_mat.values),
                     index=rho_mat.index, columns=rho_mat.columns)
fig, ax = plt.subplots(figsize=(6, 5.5))
sns.heatmap(rho_mat, cmap="RdBu_r", center=0, annot=annot, fmt="",
            cbar_kws={"label": "Spearman ρ"}, ax=ax, vmin=-0.6, vmax=0.6)
ax.set_title(f"TCGA-LUAD: TAN/NAN sig × MP score Spearman (n={len(combined)})")
fig.tight_layout()
fig.savefig(FIG / "step28_tcga_correlation_heatmap.pdf", bbox_inches="tight")
fig.savefig(FIG / "step28_tcga_correlation_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ----- 4. survival — build clinical + scores -----
print("\n[4] build survival table")
df = clin_pt.set_index("sample_barcode").join(combined, how="inner").reset_index()
df["event"] = (df["vital_status"].str.strip().str.lower() == "dead").astype(int)
df["time"] = np.where(df["event"] == 1, df["days_to_death"], df["days_to_last_follow_up"])
df = df[df["time"].notna() & (df["time"] > 0)].copy()
df["age"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce") / 365.25
df["stage_simple"] = df["ajcc_stage"].fillna("Unknown").str.extract(r"(Stage [IV]+)", expand=False).fillna("Unknown")
df["stage_num"] = df["stage_simple"].map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
df["gender_bin"] = (df["gender"].str.lower() == "male").astype(int)
print(f"  survival n={len(df)}; events={df['event'].sum()}")

# ----- 5. 4-group stratified KM -----
print("\n[5] 4-group KM: TAN-1 × MP1, TAN-4 × MP2")
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

def km4(df, sig, mp, outpath, title):
    s_med = df[sig].median(); m_med = df[mp].median()
    df = df.copy()
    df["sig_grp"] = np.where(df[sig] >= s_med, "Hi", "Lo")
    df["mp_grp"]  = np.where(df[mp] >= m_med, "Hi", "Lo")
    df["combo"] = df["sig_grp"] + "_" + df["mp_grp"]   # HiHi HiLo LoHi LoLo
    counts = df["combo"].value_counts()
    print(f"  {sig} × {mp} groups:\n{counts}")
    res = multivariate_logrank_test(df["time"], df["combo"], df["event"])
    print(f"  multivariate logrank p={res.p_value:.3g}")

    fig, ax = plt.subplots(figsize=(7, 5.5))
    palette = {"Hi_Hi": "#d62728", "Hi_Lo": "#ff7f0e", "Lo_Hi": "#2ca02c", "Lo_Lo": "#1f77b4"}
    for grp in ["Hi_Hi", "Hi_Lo", "Lo_Hi", "Lo_Lo"]:
        sub = df[df["combo"] == grp]
        if len(sub) < 5: continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub["time"], sub["event"], label=f"{grp} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=palette[grp])
    ax.set_title(f"{title}\nmultivariate logrank p={res.p_value:.3g}")
    ax.set_xlabel("days"); ax.set_ylabel("OS prob")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    fig.savefig(str(outpath).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {"sig": sig, "mp": mp, "logrank_p": res.p_value, "n_total": len(df), **counts.to_dict()}

km_rows = []
km_rows.append(km4(df, "Neu_Inflammatory", "MP1", FIG / "step28_km_neuInflammatory_mp1_4group.pdf",
                   "TCGA-LUAD: Neu_Inflammatory × MP1 4-group OS"))
km_rows.append(km4(df, "Neu_Metastatic", "MP2", FIG / "step28_km_neuMetastatic_mp2_4group.pdf",
                   "TCGA-LUAD: Neu_Metastatic × MP2 4-group OS"))
km_rows.append(km4(df, "Neu_OSM_priming", "MP1", FIG / "step28_km_neuOSMpriming_mp1_4group.pdf",
                   "TCGA-LUAD: Neu_OSM_priming × MP1 4-group OS"))
km_rows.append(km4(df, "Neu_ECM_remodeling", "MP1", FIG / "step28_km_neuECM_mp1_4group.pdf",
                   "TCGA-LUAD: Neu_ECM_remodeling × MP1 4-group OS"))
pd.DataFrame(km_rows).to_csv(RES / "step28_km_4group_logrank.csv", index=False)

# Univariate KM by single signature (median-split) — Fig 5D-equivalent forest data
print("\n[5b] univariate Cox per signature → forest plot (Fig 5D analog)")
cox_rows = []
for sg in sig_cols + mp_cols:
    sub = df[["time", "event", sg]].dropna()
    if len(sub) < 30: continue
    cph = CoxPHFitter()
    try:
        cph.fit(sub, duration_col="time", event_col="event")
        s = cph.summary.loc[sg]
        cox_rows.append({"signature": sg, "HR": s["exp(coef)"],
                         "HR_lo": s["exp(coef) lower 95%"], "HR_hi": s["exp(coef) upper 95%"],
                         "p": s["p"], "n": len(sub)})
    except Exception as e:
        print(f"  univariate fail {sg}: {e}")
cox_uni = pd.DataFrame(cox_rows).sort_values("HR", ascending=False)
cox_uni.to_csv(RES / "step28_cox_univariate.csv", index=False)
print(cox_uni.to_string(index=False))

# forest plot
fig, ax = plt.subplots(figsize=(6, 0.4 * len(cox_uni) + 1.5))
y = np.arange(len(cox_uni))[::-1]
colors = ["#d62728" if (r["HR"] > 1 and r["p"] < 0.05)
          else "#2ca02c" if (r["HR"] < 1 and r["p"] < 0.05)
          else "#888" for _, r in cox_uni.iterrows()]
ax.errorbar(cox_uni["HR"], y, xerr=[cox_uni["HR"] - cox_uni["HR_lo"], cox_uni["HR_hi"] - cox_uni["HR"]],
            fmt="o", color="black", ecolor="gray", elinewidth=1.0, capsize=3, markersize=5,
            mfc="black")
for i, (_, r) in enumerate(cox_uni.iterrows()):
    ax.scatter(r["HR"], y[i], color=colors[i], s=80, zorder=3,
               edgecolor="black", linewidth=0.7)
ax.axvline(1.0, color="black", ls="--", lw=0.5)
ax.set_yticks(y); ax.set_yticklabels(cox_uni["signature"])
ax.set_xlabel("Hazard ratio (95% CI)")
ax.set_title(f"TCGA-LUAD univariate Cox, OS (n≈{cox_uni['n'].max()})")
ax.set_xscale("log")
fig.tight_layout()
fig.savefig(FIG / "step28_cox_univariate_forest.pdf", bbox_inches="tight")
fig.savefig(FIG / "step28_cox_univariate_forest.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ----- 6. Multivariate Cox: TAN-1 + MP1 + age + stage + gender -----
print("\n[6] multivariate Cox")
def fit_multi(cols, label):
    sub = df[["time", "event"] + cols].dropna()
    if len(sub) < 50:
        print(f"  multi {label}: insufficient data ({len(sub)})"); return None
    cph = CoxPHFitter()
    cph.fit(sub, duration_col="time", event_col="event")
    out = cph.summary[["coef", "exp(coef)", "exp(coef) lower 95%",
                       "exp(coef) upper 95%", "p"]].copy()
    out["model"] = label; out["n"] = len(sub)
    out = out.reset_index().rename(columns={"covariate": "variable",
                                            "exp(coef)": "HR",
                                            "exp(coef) lower 95%": "HR_lo",
                                            "exp(coef) upper 95%": "HR_hi"})
    return out

multi_outs = []
for label, cols in [
    ("NeuInflam+MP1+covar",  ["Neu_Inflammatory", "MP1", "age", "stage_num", "gender_bin"]),
    ("NeuMet+MP2+covar",     ["Neu_Metastatic",  "MP2", "age", "stage_num", "gender_bin"]),
    ("NeuOSMpriming+MP1+covar", ["Neu_OSM_priming", "MP1", "age", "stage_num", "gender_bin"]),
    ("NeuECM+MP1+covar",     ["Neu_ECM_remodeling", "MP1", "age", "stage_num", "gender_bin"]),
    ("NeuInflam_only_covar", ["Neu_Inflammatory", "age", "stage_num", "gender_bin"]),
]:
    o = fit_multi(cols, label)
    if o is not None:
        multi_outs.append(o)
multi = pd.concat(multi_outs, ignore_index=True)
multi.to_csv(RES / "step28_cox_multivariate.csv", index=False)
print(multi.round(3).to_string(index=False))

# multi forest grouped by model
fig, axes = plt.subplots(1, len(multi_outs), figsize=(5 * len(multi_outs), 4),
                         sharex=True)
if len(multi_outs) == 1: axes = [axes]
for ax, sub in zip(axes, multi_outs):
    sub2 = sub.copy()
    y = np.arange(len(sub2))[::-1]
    colors = ["#d62728" if (r["HR"] > 1 and r["p"] < 0.05)
              else "#2ca02c" if (r["HR"] < 1 and r["p"] < 0.05)
              else "#888" for _, r in sub2.iterrows()]
    ax.errorbar(sub2["HR"], y, xerr=[sub2["HR"] - sub2["HR_lo"], sub2["HR_hi"] - sub2["HR"]],
                fmt="o", color="black", ecolor="gray", elinewidth=1.0, capsize=3,
                markersize=4)
    for i, (_, r) in enumerate(sub2.iterrows()):
        ax.scatter(r["HR"], y[i], color=colors[i], s=70, zorder=3,
                   edgecolor="black", linewidth=0.6)
    ax.axvline(1.0, color="black", ls="--", lw=0.5)
    ax.set_yticks(y); ax.set_yticklabels(sub2["variable"])
    ax.set_xscale("log")
    ax.set_xlabel("HR (95% CI)")
    ax.set_title(f"{sub['model'].iloc[0]} (n={sub['n'].iloc[0]})")
fig.suptitle("Multivariate Cox forest", y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(FIG / "step28_cox_multivariate_forest.pdf", bbox_inches="tight")
fig.savefig(FIG / "step28_cox_multivariate_forest.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"\nelapsed: {(time.time()-t0)/60:.1f} min")
print("DONE.")
