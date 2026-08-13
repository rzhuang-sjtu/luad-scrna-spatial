"""Step 13b: GSE68465 (NCI Director's Challenge LUAD, n~462) external validation.

Pipeline:
  1. Parse GSE68465 series matrix (HG-U133A / GPL96): probe × sample expression.
  2. Probe → gene symbol via GPL96 (fetched once via GEOparse).
  3. Collapse duplicates (max per gene).
  4. Parse clinical from Sample_characteristics_ch1; merge OS time + event.
  5. ssGSEA MP1-4 + EMT + Neutrophil + NETs.
  6. Univariate continuous Cox per MP; median-split KM.
  7. MP3 ↔ Neutrophil_core Spearman.

Outputs → ${WORK_ROOT}/luad_figures/fig_validation/:
  gse68465_mp_scores.csv  gse68465_km_data.csv  gse68465_cox_results.csv
  gse68465_mp3_neutrophil.csv  gse68465_summary.md
"""
from __future__ import annotations
import os, gzip, time, re
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SM = Path("${DATA_ROOT}/GSE68465/GSE68465_series_matrix.txt.gz")
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


def parse_series_matrix(path: Path):
    """Returns (expr_df: probes×samples log2 expression, sample_meta: dict[gsm]→{key:val})."""
    log("parsing series matrix")
    sample_titles = []
    sample_geo = []
    char_lines = []
    expr_lines = []
    in_table = False
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True; continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False; continue
            if in_table:
                expr_lines.append(line.rstrip("\n"))
            else:
                if line.startswith("!Sample_title\t"):
                    sample_titles = [x.strip().strip('"') for x in
                                      line.rstrip().split("\t")[1:]]
                elif line.startswith("!Sample_geo_accession\t"):
                    sample_geo = [x.strip().strip('"') for x in
                                   line.rstrip().split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1\t"):
                    char_lines.append(line.rstrip().split("\t")[1:])

    log(f"  samples: titles={len(sample_titles)} geo={len(sample_geo)} char_lines={len(char_lines)}")
    log(f"  expr rows: {len(expr_lines)}")

    # Build expression DataFrame
    header = expr_lines[0].split("\t")
    rows = []
    for ln in expr_lines[1:]:
        if ln.strip():
            rows.append(ln.split("\t"))
    expr = pd.DataFrame(rows, columns=header)
    expr = expr.rename(columns={expr.columns[0]: "ID_REF"})
    expr["ID_REF"] = expr["ID_REF"].astype(str).str.strip().str.strip('"')
    # numeric conversion
    for c in expr.columns[1:]:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")
    expr = expr.set_index("ID_REF")
    # Strip quotes from sample column names
    expr.columns = [c.strip().strip('"') for c in expr.columns]
    log(f"  expression matrix: {expr.shape}")

    # Build sample meta: dict[gsm] → {"key1": "val1", ...}
    sample_meta = {gsm: {} for gsm in sample_geo}
    for line_vals in char_lines:
        for i, val in enumerate(line_vals):
            if i >= len(sample_geo): break
            v = val.strip().strip('"')
            if ":" in v:
                k, vv = v.split(":", 1)
                sample_meta[sample_geo[i]][k.strip()] = vv.strip()
    for i, gsm in enumerate(sample_geo):
        sample_meta[gsm]["title"] = sample_titles[i] if i < len(sample_titles) else ""
    return expr, sample_meta


def main():
    t0 = time.time()
    expr, sample_meta = parse_series_matrix(SM)
    log(f"first expr cell-row: ID={expr.index[0]}  n_samples_with_value={expr.iloc[0].notna().sum()}")

    # Get GPL96 probe → gene symbol mapping via GEOparse
    log("loading GPL96 platform annotation")
    import GEOparse
    gpl = GEOparse.get_GEO("GPL96", destdir="${DATA_ROOT}/GSE68465/", silent=True)
    log(f"  GPL96 annotation: {gpl.table.shape}")
    log(f"  cols: {list(gpl.table.columns)[:8]}")
    # GPL96 has 'Gene Symbol' column
    sym_col = next((c for c in gpl.table.columns if c.lower() == "gene symbol"), None)
    log(f"  symbol column: {sym_col}")
    probe_to_gene = (gpl.table[["ID", sym_col]]
                     .rename(columns={"ID": "ID_REF", sym_col: "gene"}))
    probe_to_gene["gene"] = probe_to_gene["gene"].astype(str).str.split("///").str[0].str.strip()
    probe_to_gene = probe_to_gene[probe_to_gene["gene"].notna() & (probe_to_gene["gene"] != "")]

    # Map and collapse to gene-level (max per gene)
    log("collapsing to gene-level (max per probe per gene)")
    expr2 = expr.join(probe_to_gene.set_index("ID_REF"), how="inner")
    expr2 = expr2.dropna(subset=["gene"])
    expr_gene = expr2.groupby("gene").max().drop(columns=[]) if False else \
                expr2.groupby("gene")[expr.columns.tolist()].max()
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
        log(f"  {k}: {len(v)} genes present")

    log("running ssGSEA")
    import gseapy as gp
    ss = gp.ssgsea(
        data=expr_gene, gene_sets=gene_sets, outdir=None,
        sample_norm_method="rank", no_plot=True,
        min_size=3, max_size=10000, permutation_num=0, seed=0, threads=4,
    )
    scores = ss.res2d.pivot_table(index="Name", columns="Term",
                                   values="NES").astype(float)
    scores.index.name = "GSM"
    log(f"  scores: {scores.shape}")

    # Build clinical DataFrame from sample_meta
    log("building clinical table")
    clin = pd.DataFrame.from_dict(sample_meta, orient="index")
    log(f"  clinical fields: {[c for c in clin.columns if clin[c].notna().sum()>200][:25]}")
    log(f"  clinical shape: {clin.shape}")

    # GSE68465 clinical mapping (DC LungStudy)
    # vital_status / months_to_last_clinical_assessment / months_to_last_contact_or_death / disease_state
    vital_col = next((c for c in clin.columns
                       if "vital_status" in c.lower().replace(" ","_") or
                          "vital status" in c.lower()), None)
    log(f"  vital_status col: {vital_col}")
    # Show sample meta keys most frequent
    keycount = clin.notna().sum().sort_values(ascending=False)
    log("  top clinical fields:")
    log(keycount.head(20).to_string())

    # Survival columns (try multiple naming conventions)
    candidates = {
        "vital": ["vital_status", "vital status"],
        "months_last": ["months_to_last_clinical_assessment",
                        "months to last clinical assessment"],
        "months_death": ["months_to_last_contact_or_death",
                         "months to last contact or death"],
    }
    cols = {}
    for k, names in candidates.items():
        for n in names:
            for c in clin.columns:
                if c.lower().strip() == n.lower():
                    cols[k] = c; break
            if k in cols: break
    log(f"  survival columns: {cols}")

    # Build survival
    df = scores.join(clin, how="inner")
    if "vital" in cols:
        df["event"] = (df[cols["vital"]].astype(str).str.lower() == "dead").astype(int)
    elif "vital_status" in clin.columns:
        df["event"] = (df["vital_status"].astype(str).str.lower() == "dead").astype(int)
    # time
    if "months_last" in cols and "months_death" in cols:
        df["time_months"] = pd.to_numeric(df[cols["months_death"]], errors="coerce")
        # fallback: use months_last when months_death is NaN
        df.loc[df["time_months"].isna(), "time_months"] = pd.to_numeric(
            df.loc[df["time_months"].isna(), cols["months_last"]], errors="coerce")
    elif "months_death" in cols:
        df["time_months"] = pd.to_numeric(df[cols["months_death"]], errors="coerce")
    elif "months_last" in cols:
        df["time_months"] = pd.to_numeric(df[cols["months_last"]], errors="coerce")
    else:
        # last resort: scan any 'months' or 'days' col
        for c in clin.columns:
            if "month" in c.lower() and df.get("time_months") is None:
                df["time_months"] = pd.to_numeric(df[c], errors="coerce")
                cols["fallback_time"] = c; break

    df["time"] = df["time_months"] * 30.44
    df = df[df["time"].notna() & (df["time"] > 0) & df["event"].notna()].copy()
    log(f"  survival-ready: n={len(df)} events={int(df['event'].sum())}")

    # ----- ssGSEA results table -----
    df_out = df.reset_index()
    df_out.to_csv(OUT/"gse68465_mp_scores.csv", index=False)

    # ----- Univariate Cox per MP + median-split KM -----
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test

    cox_rows = []; km_rows = []
    for mp_col in ["MP1","MP2","MP3","MP4","Neutrophil_core","NETs_composite",
                    "EMT_Hallmark"]:
        sub = df[["time","event", mp_col]].dropna()
        if len(sub) < 30: continue
        cph = CoxPHFitter()
        cph.fit(sub, duration_col="time", event_col="event")
        s = cph.summary.loc[mp_col]
        cox_rows.append({
            "score": mp_col, "n": len(sub), "events": int(sub["event"].sum()),
            "coef": s["coef"], "HR": s["exp(coef)"],
            "HR_lo": s["exp(coef) lower 95%"], "HR_hi": s["exp(coef) upper 95%"],
            "p": s["p"], "concordance": cph.concordance_index_,
        })
        med = sub[mp_col].median()
        sub_t = sub.assign(_grp=np.where(sub[mp_col] >= med, "High", "Low"))
        lr = logrank_test(sub_t[sub_t["_grp"]=="High"]["time"],
                          sub_t[sub_t["_grp"]=="Low"]["time"],
                          sub_t[sub_t["_grp"]=="High"]["event"],
                          sub_t[sub_t["_grp"]=="Low"]["event"])
        for grp_label, sub_g in sub_t.groupby("_grp"):
            kmf = KaplanMeierFitter()
            kmf.fit(sub_g["time"], sub_g["event"], label=grp_label)
            tab = kmf.survival_function_.reset_index()
            tab.columns = ["time","surv_prob"]
            tab["score"] = mp_col; tab["group"] = grp_label
            tab["n_at_risk"] = [(sub_g["time"] >= t).sum() for t in tab["time"]]
            tab["median_cut"] = float(med)
            tab["logrank_p"] = float(lr.p_value)
            km_rows.append(tab)

    cox_df = pd.DataFrame(cox_rows)
    cox_df.to_csv(OUT/"gse68465_cox_results.csv", index=False)
    pd.concat(km_rows, ignore_index=True).to_csv(OUT/"gse68465_km_data.csv", index=False)
    log("Cox results:")
    log(cox_df.round(4).to_string(index=False))

    # MP3 ↔ Neutrophil
    sub = df[["MP3","Neutrophil_core","NETs_composite","EMT_Hallmark"]].dropna()
    rho, p = spearmanr(sub["MP3"], sub["Neutrophil_core"])
    rho_emt, p_emt = spearmanr(sub["MP3"], sub["EMT_Hallmark"])
    rho_nets, p_nets = spearmanr(sub["MP3"], sub["NETs_composite"])
    pd.DataFrame({
        "pair": ["MP3 vs Neutrophil_core", "MP3 vs EMT_Hallmark", "MP3 vs NETs_composite"],
        "rho":  [rho, rho_emt, rho_nets],
        "p":    [p, p_emt, p_nets],
        "n":    [len(sub)] * 3,
    }).to_csv(OUT/"gse68465_mp3_neutrophil.csv", index=False)
    log(f"MP3 ↔ Neut rho={rho:.3f} p={p:.2e} (n={len(sub)})")
    log(f"MP3 ↔ EMT  rho={rho_emt:.3f} p={p_emt:.2e}")
    log(f"MP3 ↔ NETs rho={rho_nets:.3f} p={p_nets:.2e}")

    # Summary
    with open(OUT/"gse68465_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 13b — GSE68465 LUAD survival validation\n\n")
        f.write(f"- Samples after probe→gene + clinical merge: n={len(df)}\n")
        f.write(f"- Events (deaths): {int(df['event'].sum())}\n\n")
        f.write("## Univariate Cox per score\n\n")
        f.write(cox_df.round(4).to_markdown(index=False) + "\n\n")
        f.write(f"## Cross-validation correlations\n\n"
                f"- MP3 ↔ Neutrophil_core: rho={rho:.3f}, p={p:.2e}\n"
                f"- MP3 ↔ EMT_Hallmark:    rho={rho_emt:.3f}, p={p_emt:.2e}\n"
                f"- MP3 ↔ NETs_composite:  rho={rho_nets:.3f}, p={p_nets:.2e}\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
