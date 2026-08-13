#!/usr/bin/env python
"""
H6 — recompute four analyses of 7p11.2 / EGFR co-amplification confounding (7p11.2 / EGFR co-amplification confounding).

Archive 1.4–1.6 holds these numbers, but intermediate tables lived in a temporary directory that has since been cleared.
Figures must be file-driven; recompute and write to disk here.

Four analyses:
  (1) SEC61G expression vs SEC61G/EGFR GISTIC copy number
  (2) OS Cox adjusting for 7p11.2 copy number + age/stage/sex
  (3) DepMap Chronos vs copy number and EGFR amplification status (already available under depmap_lineage/)
  (4) CopyKAT 7p status (already available; see archive 1.1)

Key quantities: copy-number correlation 0.984 but expression correlation only 0.178 — co-amplification does not yield co-expression,
which is the decisive evidence separating SEC61G from a mere passenger of 7p gain.

Output: results/cn_confounding/
"""
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

TPM = "${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv"
CLIN = "${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv"
CN = "${WORK_ROOT}/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"
CNT = "${WORK_ROOT}/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
OUT = "${PROJECT_ROOT}/results/cn_confounding"

# Genes on 7p11.2 co-amplified with SEC61G, plus controls elsewhere on 7p
LOCUS = ["SEC61G", "EGFR", "LANCL2", "VOPP1", "SEC61A1", "SEC61B", "SRSF9",
         "ANGPTL4"]


def main():
    os.makedirs(OUT, exist_ok=True)

    print("Loading GISTIC continuous copy number ...", flush=True)
    cn = pd.read_csv(CN, sep="\t", index_col=0)
    cn.index = cn.index.astype(str).str.upper()
    cn = cn[~cn.index.duplicated()]
    print(f"  {cn.shape}", flush=True)

    print("Loading TPM ...", flush=True)
    tpm = pd.read_csv(TPM, index_col=0)
    tpm.index = tpm.index.astype(str).str.upper()
    tpm = tpm[~tpm.index.duplicated()]
    print(f"  {tpm.shape}", flush=True)

    # Align sample barcodes to the first 15 characters (TCGA-XX-XXXX-01)
    def key(cols):
        return pd.Index([c[:15].replace(".", "-") for c in cols])
    cn.columns = key(cn.columns)
    tpm.columns = key(tpm.columns)
    common = cn.columns.intersection(tpm.columns)
    common = common[~common.duplicated()]
    print(f"shared samples {len(common)}", flush=True)
    common = pd.Index(common).drop_duplicates()
    cn = cn.loc[:, ~cn.columns.duplicated()][common]
    tpm = tpm.loc[:, ~tpm.columns.duplicated()][common]

    genes = [g for g in LOCUS if g in cn.index and g in tpm.index]
    print(f"available genes {genes}", flush=True)

    rows = []
    for g in genes:
        if g == "SEC61G":
            continue
        rows.append(dict(
            gene=g,
            cn_vs_SEC61G_cn=cn.loc["SEC61G"].corr(cn.loc[g], method="spearman"),
            expr_vs_SEC61G_expr=np.log1p(tpm.loc["SEC61G"]).corr(
                np.log1p(tpm.loc[g]), method="spearman"),
            own_expr_vs_own_cn=np.log1p(tpm.loc[g]).corr(cn.loc[g],
                                                         method="spearman")))
    C = pd.DataFrame(rows)
    C.loc[len(C)] = dict(gene="SEC61G", cn_vs_SEC61G_cn=1.0,
                         expr_vs_SEC61G_expr=1.0,
                         own_expr_vs_own_cn=np.log1p(tpm.loc["SEC61G"]).corr(
                             cn.loc["SEC61G"], method="spearman"))
    C.to_csv(f"{OUT}/cn_expr_correlations.csv", index=False)
    print("\n[Copy-number correlation vs expression correlation] — whether co-amplification yields co-expression")
    print(C.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))

    clin = pd.read_csv(CLIN)
    clin["key"] = clin.sample_barcode.astype(str).str[:15]
    clin = clin.drop_duplicates("key").set_index("key")
    # Truncating to 15 characters can create duplicate barcodes (multiple aliquots per patient); must deduplicate,
    # else reindex raises duplicate labels and the same person enters the risk set more than once
    idx = pd.Index(clin.index.intersection(common)).drop_duplicates()
    cl = clin.loc[idx]

    d = pd.DataFrame(index=idx)
    d["event"] = (cl.vital_status.astype(str).str.lower() == "dead").astype(int)
    d["time"] = np.where(d.event == 1,
                         pd.to_numeric(cl.days_to_death, errors="coerce"),
                         pd.to_numeric(cl.days_to_last_follow_up, errors="coerce"))
    age = pd.to_numeric(cl.age_at_diagnosis, errors="coerce")
    d["age"] = age / 365.25 if age.median() > 200 else age
    d["male"] = (cl.gender.astype(str).str.lower() == "male").astype(int)
    st = cl.ajcc_stage.astype(str).str.upper().str.replace("STAGE ", "", regex=False)
    d["stage_late"] = st.str.replace(r"[AB]$", "", regex=True).str.startswith(
        ("III", "IV")).astype(int)
    for g in genes:
        v = np.log1p(tpm.loc[g, idx].astype(float))
        d[f"e_{g}"] = (v - v.mean()) / v.std()
        c = cn.loc[g, idx].astype(float)
        d[f"c_{g}"] = (c - c.mean()) / c.std()
    d = d[d.time.notna() & (d.time > 0)].dropna(subset=["age"])
    print(f"\n{len(d)} analysable cases, {int(d.event.sum())} events", flush=True)

    def cox(cols):
        s = d[["time", "event"] + cols].dropna()
        c = CoxPHFitter().fit(s, "time", "event")
        return c, s

    specs = [
        ("SEC61G expression alone", ["e_SEC61G"]),
        ("+ age/stage/sex", ["e_SEC61G", "age", "stage_late", "male"]),
        ("+ SEC61G copy number", ["e_SEC61G", "c_SEC61G"]),
        ("+ 7p11.2 CN (EGFR locus)", ["e_SEC61G", "c_EGFR"]),
        ("+ 7p11.2 CN + clinical",
         ["e_SEC61G", "c_EGFR", "age", "stage_late", "male"]),
    ]
    rows = []
    for name, cols in specs:
        c, s = cox(cols)
        ci = c.confidence_intervals_.loc["e_SEC61G"]
        rows.append(dict(model=name, HR=float(np.exp(c.params_["e_SEC61G"])),
                         CI_low=float(np.exp(ci.iloc[0])),
                         CI_high=float(np.exp(ci.iloc[1])),
                         p=float(c.summary.loc["e_SEC61G", "p"]),
                         n=len(s), events=int(s.event.sum())))
    M = pd.DataFrame(rows)
    M.to_csv(f"{OUT}/sec61g_cox_ladder.csv", index=False)
    print("\n[OS Cox for SEC61G expression, stepwise adding copy number and clinical covariates]")
    print(M.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    rows = []
    for g in genes:
        cols = [f"e_{g}", "c_EGFR", "age", "stage_late", "male"]
        s = d[["time", "event"] + cols].dropna()
        c = CoxPHFitter().fit(s, "time", "event")
        ci = c.confidence_intervals_.loc[f"e_{g}"]
        rows.append(dict(gene=g, HR=float(np.exp(c.params_[f"e_{g}"])),
                         CI_low=float(np.exp(ci.iloc[0])),
                         CI_high=float(np.exp(ci.iloc[1])),
                         p=float(c.summary.loc[f"e_{g}", "p"]), n=len(s)))
    L = pd.DataFrame(rows).sort_values("p")
    L["bonferroni_pass"] = L.p < 0.05 / len(L)
    L.to_csv(f"{OUT}/locus_specificity.csv", index=False)
    print("\n[Locus specificity: same model (+7p11.2 CN + clinical), swapping genes one at a time]")
    print(L.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(f"\nWriting {OUT}/", flush=True)


if __name__ == "__main__":
    main()
