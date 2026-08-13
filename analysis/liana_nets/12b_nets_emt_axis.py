"""Step 12b: NETs-EMT axis — computational validation.

Three parts:
  B. TCGA NETs ssGSEA → correlate with MP3/EMT/Neutrophil + joint survival
  C. GSE253013 Granulocyte analysis (author-annotated 4,933 cells)
  D. NETs → STING/cGAS → IFN-I mechanism check in TCGA

Skips Step A (scRNA NETs scoring) because only 2/8 core NET genes survive
the 9881-HVG filter; bulk TCGA has full coverage.
"""
from __future__ import annotations
import os, gc, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

TPM = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv")
CLIN = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
MP_SSGSEA = Path("${WORK_ROOT}/luad_figures/fig3/tcga_luad_mp_ssgsea.csv.gz")
NEUT_SCORES = Path("${WORK_ROOT}/luad_figures/fig3/tcga_neutrophil_scores.csv")
GMT = Path.home() / "luad/data/gmt/MSigDB_Hallmark_2020.gmt"
COPYKAT = Path.home() / "luad/data/processed/luad_copykat.h5ad"
MP_SCORES = Path.home() / "luad/results/step7_mp_cell_scores.csv"

RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig_nets")
FIG.mkdir(parents=True, exist_ok=True)

# -- NETs panel (Yang Nat Commun 2023 + Jiang Front Immunol 2022 consensus) --
NETS_PANELS = {
    "NETs_core":     ["PADI4", "MPO", "ELANE", "CTSG", "PRTN3", "DEFA1", "DEFA3"],
    "NETs_release":  ["HMGB1", "H3F3A", "DNASE1L3"],
    "NADPH_oxidase": ["CYBB", "NCF1", "NCF2", "NCF4"],
    "NETs_receptor": ["TLR4", "TLR9", "RIPK3", "PAD4"],
    "NETs_markers":  ["S100A8", "S100A9", "S100A12", "CAMP", "LCN2", "MMP9", "OLFM4"],
}
NETS_COMPOSITE = (NETS_PANELS["NETs_core"] + NETS_PANELS["NETs_release"]
                   + NETS_PANELS["NADPH_oxidase"] + NETS_PANELS["NETs_markers"])

STING_CGAS = ["TMEM173", "STING1", "MB21D1", "CGAS"]
IFN_I_SIG  = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
              "OAS1", "OAS2", "OAS3", "IFI6", "STAT1", "IRF7"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def hallmark_emt() -> list[str]:
    with open(GMT) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0].lower().startswith("epithelial mesenchymal"):
                return parts[2:]
    return []


def mean_z(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean of per-gene z-scored log-expression (sample axis)."""
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns, name="score")
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0)


def main() -> None:
    t0 = time.time()

    # B. TCGA NETs ssGSEA + correlations + survival
    log("loading TCGA TPM")
    tpm = pd.read_csv(TPM, index_col=0)
    if not tpm.index.is_unique:
        tpm = tpm.groupby(tpm.index).max()
    log(f"  TPM: {tpm.shape}")
    clin = pd.read_csv(CLIN)
    clin_pt = clin[clin["sample_type"] == "Primary Tumor"].copy()
    pt_samples = [s for s in tpm.columns if s in set(clin_pt["sample_barcode"])]
    tpm = tpm[pt_samples]
    log(f"  primary tumor samples: {len(pt_samples)}")

    expr = np.log2(tpm + 1.0).astype("float32")

    # --- ssGSEA for NETs composite + STING/IFN ---
    log("ssGSEA: NETs composite, STING/cGAS, IFN-I")
    import gseapy as gp
    gs = {
        "NETs_composite": [g for g in NETS_COMPOSITE if g in expr.index],
        "NETs_core":      [g for g in NETS_PANELS["NETs_core"] if g in expr.index],
        "STING_cGAS":     [g for g in STING_CGAS if g in expr.index],
        "IFN_I":          [g for g in IFN_I_SIG if g in expr.index],
        "EMT_Hallmark":   [g for g in hallmark_emt() if g in expr.index],
    }
    for k, v in gs.items():
        log(f"  {k}: {len(v)} genes present")

    ss = gp.ssgsea(
        data=expr, gene_sets=gs, outdir=None,
        sample_norm_method="rank", no_plot=True,
        min_size=3, max_size=5000, permutation_num=0, seed=0, threads=8,
    )
    scores = ss.res2d.pivot_table(index="Name", columns="Term", values="NES").astype(float)
    scores.index.name = "sample_barcode"

    # Also mean_z variants for robustness
    scores["NETs_meanZ"] = mean_z(expr, gs["NETs_composite"])
    scores["STING_meanZ"] = mean_z(expr, gs["STING_cGAS"])
    scores["IFN_meanZ"] = mean_z(expr, gs["IFN_I"])

    # Load MP ssGSEA + Neutrophil score
    mp = pd.read_csv(MP_SSGSEA, index_col=0)
    mp.columns = [f"MP_{c}" for c in mp.columns]
    neut = pd.read_csv(NEUT_SCORES, index_col=0)

    # Drop conflicting cols from older Neutrophil CSV to avoid join collision
    neut = neut[[c for c in ["Neutrophil_core", "TAN_proTumoral"] if c in neut.columns]]
    merged = scores.join([mp, neut], how="inner")
    log(f"  merged: {merged.shape}")
    merged.to_csv(FIG / "tcga_nets_scores_full.csv.gz", compression="gzip")

    # --- Spearman correlations ---
    log("Spearman correlations")
    corr_rows = []
    for a_col in ["NETs_composite", "NETs_core", "NETs_meanZ"]:
        for b_col in ["MP_MP1", "MP_MP2", "MP_MP3", "MP_MP4",
                       "EMT_Hallmark", "Neutrophil_core", "TAN_proTumoral",
                       "STING_cGAS", "STING_meanZ", "IFN_I", "IFN_meanZ"]:
            if b_col not in merged.columns: continue
            sub = merged[[a_col, b_col]].dropna()
            if len(sub) < 5: continue
            rho, p = spearmanr(sub[a_col], sub[b_col])
            corr_rows.append({"score_A": a_col, "score_B": b_col,
                              "spearman_rho": rho, "p": p, "n": len(sub)})
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(FIG / "nets_score_mp_correlation.csv", index=False)
    corr_df.to_csv(RES / "step12b_nets_tcga_correlations.csv", index=False)
    log("NETs correlations (top):")
    log(corr_df[corr_df["score_A"] == "NETs_composite"]
        .sort_values("spearman_rho", ascending=False).round(4).to_string(index=False))

    # EMT-specific correlation file
    emt_rows = corr_df[corr_df["score_B"] == "EMT_Hallmark"].copy()
    emt_rows.to_csv(FIG / "nets_emt_correlation.csv", index=False)

    # STING/IFN mechanism table
    sting_rows = corr_df[corr_df["score_B"].str.contains("STING|IFN")].copy()
    sting_rows.to_csv(FIG / "nets_sting_ifn_correlation.csv", index=False)

    # --- Combined survival: NETs × MP3 ---
    log("2×2 joint KM: NETs × MP3")
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test

    surv = clin_pt.set_index("sample_barcode").join(merged, how="inner").reset_index()
    surv["event"] = (surv["vital_status"].str.strip().str.lower() == "dead").astype(int)
    surv["time"] = np.where(surv["event"] == 1, surv["days_to_death"],
                              surv["days_to_last_follow_up"])
    surv = surv[surv["time"].notna() & (surv["time"] > 0)].copy()
    log(f"  surv: n={len(surv)} events={int(surv['event'].sum())}")

    nmed = surv["NETs_composite"].median()
    mmed = surv["MP_MP3"].median()
    surv["NETs_grp"] = np.where(surv["NETs_composite"] >= nmed, "NETs_hi", "NETs_lo")
    surv["MP3_grp"] = np.where(surv["MP_MP3"] >= mmed, "MP3hi", "MP3lo")
    surv["joint"] = surv["MP3_grp"] + "_" + surv["NETs_grp"]

    km_rows = []
    for g, sub in surv.groupby("joint"):
        kmf = KaplanMeierFitter()
        kmf.fit(sub["time"], sub["event"], label=g)
        tab = kmf.survival_function_.reset_index()
        tab.columns = ["time", "surv_prob"]
        tab["group"] = g
        tab["n_at_risk"] = [(sub["time"] >= t).sum() for t in tab["time"]]
        tab["n_total"] = len(sub)
        tab["events_total"] = int(sub["event"].sum())
        km_rows.append(tab)
    pd.concat(km_rows, ignore_index=True).to_csv(
        FIG / "nets_mp3_km.csv", index=False)

    # pairwise + overall logrank
    pair = []
    grps = sorted(surv["joint"].unique())
    for i in range(len(grps)):
        for j in range(i+1, len(grps)):
            a = surv[surv["joint"] == grps[i]]
            b = surv[surv["joint"] == grps[j]]
            lr = logrank_test(a["time"], b["time"], a["event"], b["event"])
            pair.append({"group_A": grps[i], "group_B": grps[j],
                         "n_A": len(a), "n_B": len(b),
                         "chi2": lr.test_statistic, "p": lr.p_value})
    overall = multivariate_logrank_test(surv["time"], surv["joint"], surv["event"])
    pair.append({"group_A": "ALL", "group_B": "ALL", "n_A": len(surv),
                 "n_B": len(surv), "chi2": overall.test_statistic,
                 "p": overall.p_value})
    pd.DataFrame(pair).to_csv(FIG / "nets_mp3_logrank.csv", index=False)

    # NETs-alone median-split KM
    lr_nets = logrank_test(surv[surv["NETs_grp"]=="NETs_hi"]["time"],
                            surv[surv["NETs_grp"]=="NETs_lo"]["time"],
                            surv[surv["NETs_grp"]=="NETs_hi"]["event"],
                            surv[surv["NETs_grp"]=="NETs_lo"]["event"])
    log(f"  NETs-alone log-rank chi2={lr_nets.test_statistic:.2f} p={lr_nets.p_value:.4f}")

    # --- Multivariate Cox: NETs + MP3 + stage + age + gender ---
    log("multivariate Cox")
    surv["age"] = pd.to_numeric(surv["age_at_diagnosis"], errors="coerce") / 365.25
    surv["stage_num"] = (surv["ajcc_stage"].fillna("Unknown")
                         .str.extract(r"(Stage [IV]+)", expand=False)
                         .map({"Stage I":1,"Stage II":2,"Stage III":3,"Stage IV":4}))
    surv["gender_bin"] = (surv["gender"].str.lower() == "male").astype(int)
    cox_df = surv[["time","event","age","stage_num","gender_bin",
                    "MP_MP3","NETs_composite"]].dropna().copy()
    # Z-score continuous covariates
    for c in ["MP_MP3","NETs_composite"]:
        cox_df[c] = (cox_df[c] - cox_df[c].mean()) / cox_df[c].std()
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col="time", event_col="event")
    mv = cph.summary[["coef","exp(coef)","exp(coef) lower 95%",
                       "exp(coef) upper 95%","p"]].copy()
    mv.columns = ["coef","HR","HR_lo","HR_hi","p"]
    mv["concordance"] = cph.concordance_index_
    mv.to_csv(FIG / "nets_mp3_cox.csv")
    log(mv.round(4).to_string())

    # C. GSE253013 Granulocyte analysis
    log("loading copykat h5ad for GSE253013 Granulocyte analysis")
    a = sc.read_h5ad(COPYKAT)
    gse = a[a.obs["dataset"] == "GSE253013"].copy()
    log(f"  GSE253013 cells: {gse.n_obs}")

    gr_mask = (gse.obs["celltype_original"] == "Granulocytes")
    log(f"  Granulocytes: {int(gr_mask.sum())}")

    # Top 50 expressed genes in Granulocytes (mean log1p)
    gran = gse[gr_mask].copy()
    # Use lognorm layer for proper log-norm expression
    gran.X = gran.layers["lognorm"]
    X_arr = gran.X.toarray() if hasattr(gran.X, "toarray") else np.asarray(gran.X)
    mean_expr = pd.Series(X_arr.mean(axis=0), index=gran.var_names)
    pct_expr  = pd.Series((X_arr > 0).mean(axis=0), index=gran.var_names)
    top50 = mean_expr.sort_values(ascending=False).head(50)
    top50_df = pd.DataFrame({"gene": top50.index,
                              "mean_log1p": top50.values,
                              "pct_expressing": pct_expr[top50.index].values})
    top50_df.to_csv(FIG / "gse253013_granulocyte_profile.csv", index=False)
    log("top-15 expressed genes in Granulocytes:")
    log(top50_df.head(15).round(3).to_string(index=False))

    # Check NETs-related genes that ARE in var
    nets_in_var = [g for g in NETS_COMPOSITE if g in gran.var_names]
    nets_in_mc  = [g for g in NETS_COMPOSITE
                   if g not in gran.var_names and g in gran.obsm.get("marker_counts",
                                                                       pd.DataFrame()).columns]
    log(f"  NETs genes in var: {nets_in_var}")
    log(f"  NETs genes in mc:  {nets_in_mc}")
    if nets_in_var:
        net_expr = pd.Series(
            X_arr[:, [list(gran.var_names).index(g) for g in nets_in_var]].mean(axis=0),
            index=nets_in_var, name="mean_log1p_in_Granulocyte")
        non_gr_mask = (gse.obs["celltype_original"] != "Granulocytes").values
        nongr = gse[non_gr_mask].copy()
        nongr.X = nongr.layers["lognorm"]
        X2 = nongr.X.toarray() if hasattr(nongr.X, "toarray") else np.asarray(nongr.X)
        non_gr_expr = pd.Series(
            X2[:, [list(nongr.var_names).index(g) for g in nets_in_var]].mean(axis=0),
            index=nets_in_var, name="mean_log1p_nonGranulocyte")
        nets_cmp = pd.concat([net_expr, non_gr_expr], axis=1)
        nets_cmp["fold"] = (nets_cmp["mean_log1p_in_Granulocyte"] /
                             nets_cmp["mean_log1p_nonGranulocyte"].replace(0, np.nan))
        nets_cmp.to_csv(FIG / "gse253013_nets_gene_comparison.csv")
        log("NETs gene expression in Granulocyte vs non-Granulocyte:")
        log(nets_cmp.round(3).to_string())

    # Per-patient Granulocyte % of immune cells and correlation with mean MP3
    log("per-patient Granulocyte % vs mean MP3")
    IMMUNE = {"T lymphocytes","T cells","Myeloid cells","Myeloid","B lymphocytes",
              "B cells","NK cells","MAST cells","Granulocytes","CD45+"}
    imm_mask = gse.obs["celltype_original"].isin(IMMUNE)
    imm = gse[imm_mask].copy()
    gran_pct = (imm.obs.groupby("patient_id")["celltype_original"]
                   .apply(lambda s: (s == "Granulocytes").mean())
                   .rename("Granulocyte_pct_immune"))
    imm_n = imm.obs.groupby("patient_id").size().rename("n_immune_cells")

    # MP3 per patient (from step7 malignant scores)
    mp_cell = pd.read_csv(MP_SCORES, index_col=0)
    # map patient_id from the full adata barcodes
    # step7 index is barcode; we need patient_id join via gse obs (malignant are elsewhere)
    # use copykat h5ad obs for barcode→patient
    bc_to_patient = a.obs["patient_id"].astype(str).to_dict()
    mp_cell["patient_id"] = mp_cell.index.map(bc_to_patient.get)
    # Only GSE253013 patients (intersect)
    gse_patients = set(gse.obs["patient_id"].astype(str).unique())
    mp_gse = mp_cell[mp_cell["patient_id"].isin(gse_patients)]
    mp_pat_mp3 = mp_gse.groupby("patient_id")["MP3_score"].mean().rename("mean_MP3")
    n_malignant = mp_gse.groupby("patient_id").size().rename("n_malignant")

    per_patient = pd.concat([gran_pct, imm_n, mp_pat_mp3, n_malignant], axis=1)
    per_patient.to_csv(FIG / "gse253013_patient_mp3_neutrophil.csv")
    log(per_patient.round(3).to_string())
    if len(per_patient.dropna()) >= 3:
        rho, p = spearmanr(per_patient["Granulocyte_pct_immune"].dropna(),
                            per_patient["mean_MP3"].dropna())
        log(f"Spearman(Granulocyte_pct, mean_MP3) rho={rho:.3f} p={p:.4f}")
    else:
        log("not enough patients for correlation")

    # Tissue-type Granulocyte % comparison in GSE253013
    tissue_gr = (gse.obs[gse.obs["celltype_original"].isin(IMMUNE)]
                 .assign(is_gran=lambda d: (d["celltype_original"]=="Granulocytes").astype(int))
                 .groupby("tissue_type")["is_gran"].mean().rename("Granulocyte_pct_of_immune"))
    tissue_gr.to_csv(FIG / "gse253013_granulocyte_by_tissue.csv")
    log("Granulocyte % of immune by tissue_type:")
    log(tissue_gr.round(3).to_string())

    # Save summary md
    with open(RES / "step12b_nets_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 12b — NETs-EMT axis validation\n\n")
        f.write("## A. scRNA-seq NETs gene coverage (skipped)\n\n")
        f.write("Only 2/8 NETs core genes (HMGB1, H3F3A) survive the 9881-HVG "
                "filter in luad_copykat.h5ad. The key neutrophil-specific "
                "enzymes (PADI4, MPO, ELANE, CTSG) are absent from var. "
                "Analysis moved to TCGA bulk where coverage is complete.\n\n")

        f.write("## B. TCGA NETs correlations\n\n")
        f.write(f"- n={merged.shape[0]} primary tumors with ssGSEA scores.\n")
        f.write(f"- Gene-set coverage: NETs composite "
                f"{len(gs['NETs_composite'])}, NETs core {len(gs['NETs_core'])}, "
                f"STING/cGAS {len(gs['STING_cGAS'])}, IFN-I {len(gs['IFN_I'])}.\n\n")
        f.write("### NETs composite ↔ other scores\n\n")
        nets_composite = corr_df[corr_df['score_A']=='NETs_composite'].copy()
        f.write(nets_composite.round(4).to_markdown(index=False) + "\n\n")

        f.write(f"### NETs × MP3 joint survival (n={len(surv)} events={int(surv['event'].sum())})\n\n")
        f.write(f"- NETs-alone log-rank chi2={lr_nets.test_statistic:.2f}, "
                f"p={lr_nets.p_value:.4f}\n")
        pair_df = pd.DataFrame(pair)
        f.write("\nPairwise + overall 2×2 log-rank:\n\n")
        f.write(pair_df.round(5).to_markdown(index=False) + "\n\n")
        f.write("### Multivariate Cox (z-scored NETs + MP3 + age + stage + gender)\n\n")
        f.write(mv.round(4).to_markdown() + "\n\n")

        f.write("## C. GSE253013 Granulocyte analysis\n\n")
        f.write(f"- Granulocytes labeled by original author: {int(gr_mask.sum())}\n")
        f.write(f"- Tissue distribution:\n\n")
        f.write(gse.obs[gse.obs['celltype_original']=='Granulocytes']
                ['tissue_type'].value_counts().to_frame("n_cells").to_markdown() + "\n\n")

        f.write("Top 15 expressed genes in Granulocyte (mean log1p):\n\n")
        f.write(top50_df.head(15).round(3).to_markdown(index=False) + "\n\n")

        if nets_in_var:
            f.write("NETs genes in Granulocyte vs non-Granulocyte:\n\n")
            f.write(nets_cmp.round(3).to_markdown() + "\n\n")

        f.write("Per-patient table:\n\n")
        f.write(per_patient.round(3).to_markdown() + "\n\n")

        f.write("Granulocyte % of immune cells by tissue:\n\n")
        f.write(tissue_gr.round(3).to_frame("Gran_pct_immune").to_markdown() + "\n\n")

        f.write("## D. STING/cGAS/IFN-I mechanism\n\n")
        sting_sub = corr_df[corr_df['score_B'].str.contains('STING|IFN')]
        f.write(sting_sub.round(4).to_markdown(index=False) + "\n")

    # also save a separate granulocyte patient csv
    per_patient.to_csv(RES / "step12b_gse253013_granulocyte.csv")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
