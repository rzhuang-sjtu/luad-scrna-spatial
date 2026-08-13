"""Step 13c: GSE31210 (Tokyo, n=246, stage I-II LUAD) external validation.

GPL570 (HG-U133 Plus 2.0). Already downloaded SOFT via GEOparse to
${DATA_ROOT}/GSE31210/GSE31210_family.soft.gz.

Pipeline mirrors step 13b but uses days_before_death/relapse fields.
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SOFT = Path("${DATA_ROOT}/GSE31210/GSE31210_family.soft.gz")
SIG = Path.home() / "luad/results/step6_mp_signatures_top100.csv"
GMT = Path.home() / "luad/data/gmt/MSigDB_Hallmark_2020.gmt"
OUT = Path("${WORK_ROOT}/luad_figures/fig_validation")
OUT.mkdir(parents=True, exist_ok=True)

NEUT_CORE = ["CSF3R","FCGR3B","CXCR1","CXCR2","S100A8","S100A9","MMP9","ELANE","CEACAM8"]
NETS_COMPOSITE = ["PADI4","MPO","ELANE","CTSG","PRTN3","HMGB1","H3F3A","DNASE1L3",
                  "CYBB","NCF1","NCF2","NCF4","S100A8","S100A9","S100A12","CAMP","LCN2","MMP9"]
TOPN_MP = 50


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def hallmark_emt() -> list[str]:
    with open(GMT) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0].lower().startswith("epithelial mesenchymal"):
                return parts[2:]
    return []


def main():
    t0 = time.time()
    log(f"loading GSE31210 SOFT")
    import GEOparse
    gse = GEOparse.get_GEO(filepath=str(SOFT), silent=True)
    log(f"  n_samples: {len(gse.gsms)}")

    # Build sample × probe matrix
    log("building expression matrix")
    expr_dict = {}
    for gsm_name, gsm in gse.gsms.items():
        t = gsm.table.set_index("ID_REF")["VALUE"]
        expr_dict[gsm_name] = t
    expr = pd.DataFrame(expr_dict)
    expr = expr.astype("float32")
    log(f"  probe matrix: {expr.shape}")

    # Probe → gene via GPL570 (already in gse object)
    log("probe → gene mapping")
    gpl = list(gse.gpls.values())[0]
    sym_col = next((c for c in gpl.table.columns if c.lower() == "gene symbol"), None)
    log(f"  symbol col: {sym_col}")
    p2g = gpl.table[["ID", sym_col]].rename(columns={"ID":"ID_REF", sym_col:"gene"})
    p2g["gene"] = p2g["gene"].astype(str).str.split("///").str[0].str.strip()
    p2g = p2g[p2g["gene"].notna() & (p2g["gene"] != "") & (p2g["gene"] != "nan")]
    p2g = p2g.set_index("ID_REF")

    expr2 = expr.join(p2g, how="inner")
    expr_gene = expr2.groupby("gene")[expr.columns.tolist()].max()
    log(f"  gene-level matrix: {expr_gene.shape}")

    # ssGSEA
    log("preparing gene sets")
    sig_df = pd.read_csv(SIG)
    sig_df = sig_df[(sig_df["MP"].isin(["MP1","MP2","MP3","MP4"])) &
                    (sig_df["rank"] <= TOPN_MP)]
    gene_sets = {f"MP{i}": [g for g in sig_df[sig_df["MP"]==f"MP{i}"]["gene"].tolist()
                            if g in expr_gene.index] for i in range(1,5)}
    gene_sets["EMT_Hallmark"] = [g for g in hallmark_emt() if g in expr_gene.index]
    gene_sets["Neutrophil_core"] = [g for g in NEUT_CORE if g in expr_gene.index]
    gene_sets["NETs_composite"] = [g for g in NETS_COMPOSITE if g in expr_gene.index]
    for k, v in gene_sets.items():
        log(f"  {k}: {len(v)} genes")

    log("running ssGSEA")
    import gseapy as gp
    ss = gp.ssgsea(
        data=expr_gene, gene_sets=gene_sets, outdir=None,
        sample_norm_method="rank", no_plot=True,
        min_size=3, max_size=10000, permutation_num=0, seed=0, threads=4,
    )
    scores = ss.res2d.pivot_table(index="Name", columns="Term", values="NES").astype(float)
    scores.index.name = "GSM"
    log(f"  scores: {scores.shape}")

    # Build clinical from GSM characteristics
    log("building clinical table")
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        d = {"GSM": gsm_name}
        for line in gsm.metadata.get("characteristics_ch1", []):
            if ":" in line:
                k, v = line.split(":", 1)
                d[k.strip()] = v.strip()
        rows.append(d)
    clin = pd.DataFrame(rows).set_index("GSM")
    log(f"  clinical shape: {clin.shape}")

    df = scores.join(clin, how="inner")

    # Survival fields:
    #   death (alive/dead) → event
    #   days before death/censor → time
    df["event"] = (df["death"].astype(str).str.lower() == "dead").astype(int)
    df["time"] = pd.to_numeric(df["days before death/censor"], errors="coerce")
    df["relapse_event"] = (df["relapse"].astype(str).str.lower() == "relapsed").astype(int)
    df["relapse_time"] = pd.to_numeric(df["days before relapse/censor"], errors="coerce")
    # Tumor only: keep "primary lung tumor"
    df = df[df["tissue"].astype(str).str.contains("primary", case=False, na=False)].copy()
    log(f"  primary tumor samples: {len(df)}")
    df = df[df["time"].notna() & (df["time"] > 0)].copy()
    log(f"  survival ready: n={len(df)} events={int(df['event'].sum())}")

    df.reset_index().to_csv(OUT/"gse31210_mp_scores.csv", index=False)

    # Cox + KM (OS and relapse)
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test

    def survival_block(df_in, time_col, event_col, suffix):
        sub_all = df_in[[time_col, event_col]].dropna()
        if len(sub_all) < 30 or sub_all[event_col].sum() < 5:
            log(f"  {suffix}: too few events, skipping")
            return [], []
        cox_rows = []; km_rows = []
        for mp_col in ["MP1","MP2","MP3","MP4","Neutrophil_core","NETs_composite",
                        "EMT_Hallmark"]:
            sub = df_in[[time_col, event_col, mp_col]].dropna()
            if len(sub) < 30: continue
            cph = CoxPHFitter()
            cph.fit(sub.rename(columns={time_col:"time", event_col:"event"}),
                    duration_col="time", event_col="event")
            s = cph.summary.loc[mp_col]
            cox_rows.append({
                "endpoint": suffix, "score": mp_col,
                "n": len(sub), "events": int(sub[event_col].sum()),
                "coef": s["coef"], "HR": s["exp(coef)"],
                "HR_lo": s["exp(coef) lower 95%"], "HR_hi": s["exp(coef) upper 95%"],
                "p": s["p"], "concordance": cph.concordance_index_,
            })
            med = sub[mp_col].median()
            sub_t = sub.assign(_grp=np.where(sub[mp_col]>=med, "High", "Low"))
            lr = logrank_test(sub_t[sub_t["_grp"]=="High"][time_col],
                               sub_t[sub_t["_grp"]=="Low"][time_col],
                               sub_t[sub_t["_grp"]=="High"][event_col],
                               sub_t[sub_t["_grp"]=="Low"][event_col])
            for grp_label, sub_g in sub_t.groupby("_grp"):
                kmf = KaplanMeierFitter()
                kmf.fit(sub_g[time_col], sub_g[event_col], label=grp_label)
                tab = kmf.survival_function_.reset_index()
                tab.columns = ["time","surv_prob"]
                tab["score"] = mp_col; tab["group"] = grp_label
                tab["endpoint"] = suffix
                tab["median_cut"] = float(med)
                tab["logrank_p"] = float(lr.p_value)
                km_rows.append(tab)
        return cox_rows, km_rows

    os_cox, os_km = survival_block(df, "time", "event", "OS")
    rfs_cox, rfs_km = survival_block(df, "relapse_time", "relapse_event", "RFS")

    cox_all = pd.DataFrame(os_cox + rfs_cox)
    cox_all.to_csv(OUT/"gse31210_cox_results.csv", index=False)
    pd.concat(os_km + rfs_km, ignore_index=True).to_csv(
        OUT/"gse31210_km_data.csv", index=False)
    log("Cox results:")
    log(cox_all.round(4).to_string(index=False))

    # Correlations
    sub = df[["MP3","Neutrophil_core","NETs_composite","EMT_Hallmark"]].dropna()
    rho, p = spearmanr(sub["MP3"], sub["Neutrophil_core"])
    rho_emt, p_emt = spearmanr(sub["MP3"], sub["EMT_Hallmark"])
    rho_nets, p_nets = spearmanr(sub["MP3"], sub["NETs_composite"])
    pd.DataFrame({
        "pair": ["MP3 vs Neutrophil_core","MP3 vs EMT_Hallmark","MP3 vs NETs_composite"],
        "rho": [rho, rho_emt, rho_nets],
        "p":   [p, p_emt, p_nets],
        "n":   [len(sub)]*3,
    }).to_csv(OUT/"gse31210_mp3_neutrophil.csv", index=False)
    log(f"MP3 ↔ Neut rho={rho:.3f} p={p:.2e}")
    log(f"MP3 ↔ EMT  rho={rho_emt:.3f} p={p_emt:.2e}")
    log(f"MP3 ↔ NETs rho={rho_nets:.3f} p={p_nets:.2e}")

    with open(OUT/"gse31210_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 13c — GSE31210 stage I-II LUAD validation\n\n")
        f.write(f"- Primary tumor samples: n={len(df)}\n")
        f.write(f"- OS events: {int(df['event'].sum())}; "
                f"RFS events: {int(df['relapse_event'].dropna().sum())}\n\n")
        f.write("## Cox per score (OS + RFS)\n\n")
        f.write(cox_all.round(4).to_markdown(index=False) + "\n\n")
        f.write("## Cross-validation correlations (MP3-centric)\n\n"
                f"- MP3 ↔ Neutrophil_core: rho={rho:.3f}, p={p:.2e}\n"
                f"- MP3 ↔ EMT_Hallmark:    rho={rho_emt:.3f}, p={p_emt:.2e}\n"
                f"- MP3 ↔ NETs_composite:  rho={rho_nets:.3f}, p={p_nets:.2e}\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
