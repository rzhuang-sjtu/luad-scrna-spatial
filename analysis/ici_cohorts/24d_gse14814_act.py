"""Step 24d: GSE14814 (JBR.10 trial, n=133 NSCLC) — Cox MP × ACT interaction.

Adjuvant chemotherapy (ACT, cisplatin/vinorelbine) vs OBS arms in JBR.10.
Goal: Does MP score modify ACT benefit on OS?
  - Continuous Cox: OS ~ MP_score * treatment + Stage + Histology
  - Median-split MP: HR(ACT vs OBS) within high vs low; ratio
  - ADC subset (predominant histology) + all-NSCLC (133)

Outputs (fig_treatment/):
  gse14814_mp_scores.csv          sample × MP1-4 + EMT/Neut/NETs + clinical
  gse14814_response_comparison.csv per MP: interaction p, HR(ACT|high/low), KM medians
  gse14814_boxplot_data.csv        long form for plotting
"""
from __future__ import annotations
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(os.path.expanduser("~/luad/scripts")))
from pathlib import Path
import numpy as np
import pandas as pd

from importlib import import_module
C = import_module("24_common")

GSE = "gse14814"


def main():
    C.log(f"=== {GSE.upper()} JBR.10 ACT vs OBS ===")
    import GEOparse
    gse = GEOparse.get_GEO(geo="GSE14814", destdir=str(C.EXT), silent=True,
                             annotate_gpl=True, include_data=True)
    C.log(f"n_samples: {len(gse.gsms)}")

    # 1. expr matrix probe × sample
    expr = pd.DataFrame({g: s.table.set_index("ID_REF")["VALUE"]
                         for g, s in gse.gsms.items()})
    expr = expr.astype("float32")
    C.log(f"probe matrix: {expr.shape}")

    # 2. probe → gene
    gpl = list(gse.gpls.values())[0]
    sym_col = "Gene Symbol"
    p2g = gpl.table[["ID", sym_col]].copy()
    p2g.columns = ["ID_REF", "gene"]
    p2g["gene"] = p2g["gene"].astype(str).str.split("///").str[0].str.strip()
    p2g = p2g[p2g["gene"].notna() & (p2g["gene"] != "") & (p2g["gene"] != "nan")]
    p2g = p2g.set_index("ID_REF")

    expr2 = expr.join(p2g, how="inner")
    expr_gene = expr2.groupby("gene")[expr.columns.tolist()].max()
    C.log(f"gene matrix: {expr_gene.shape}")

    # 3. ssGSEA
    gene_sets, ov = C.build_gene_sets(expr_gene.index)
    for k, n in ov.items():
        C.log(f"  {k}: {n} present")
    C.log("running ssGSEA")
    scores = C.run_ssgsea(expr_gene, gene_sets)
    C.log(f"  scores: {scores.shape}")

    # 4. clinical from characteristics
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        d = {"GSM": gsm_name}
        for c in gsm.metadata.get("characteristics_ch1", []):
            if ":" in c:
                k, v = c.split(":", 1)
                k = k.strip().lower().replace(" ", "_")
                d[k] = v.strip()
        rows.append(d)
    clin = pd.DataFrame(rows).set_index("GSM")
    C.log(f"clinical cols: {clin.columns.tolist()}")
    C.log(f"treatment counts: {clin['post_surgical_treatment'].value_counts().to_dict()}")
    C.log(f"histology: {clin['histology_type'].value_counts().to_dict()}")

    clin["os_time"] = pd.to_numeric(clin["os_time"], errors="coerce")
    clin["os_event"] = (clin["os_status"].astype(str).str.lower()
                         .map(lambda x: 1 if "dead" in x else 0)).astype(int)
    clin["treatment"] = clin["post_surgical_treatment"].map(
        {"ACT": 1, "OBS": 0}).astype(float)

    # 5. join
    score_df = scores.join(clin[["treatment", "os_time", "os_event",
                                    "stage", "histology_type", "sex", "age"]],
                              how="inner")
    score_df = score_df.dropna(subset=["treatment", "os_time", "os_event"])
    score_df.index.name = "Sample"
    score_df.to_csv(C.OUT / f"{GSE}_mp_scores.csv")
    C.log(f"after QC: {score_df.shape}")

    # 6. Cox interaction per MP
    from lifelines import CoxPHFitter
    SCORES = ["MP1", "MP2", "MP3", "MP4", "EMT_Hallmark",
               "Neutrophil_core", "NETs_composite"]
    SUBSETS = {
        "all_nsclc": score_df.index.tolist(),
        "ADC_only": score_df[score_df["histology_type"] == "ADC"].index.tolist(),
    }

    rows = []
    box_rows = []
    for sub_name, idx in SUBSETS.items():
        sub = score_df.loc[idx].copy()
        for s in SCORES:
            if s not in sub.columns:
                continue
            df = sub[[s, "treatment", "os_time", "os_event"]].dropna().copy()
            if len(df) < 20:
                continue
            df["MP_z"] = (df[s] - df[s].mean()) / df[s].std(ddof=0)
            df["MP_x_treat"] = df["MP_z"] * df["treatment"]
            try:
                cph = CoxPHFitter(penalizer=0.001)
                cph.fit(df[["MP_z", "treatment", "MP_x_treat",
                              "os_time", "os_event"]],
                        duration_col="os_time", event_col="os_event",
                        show_progress=False)
                summ = cph.summary
                interaction_p = float(summ.loc["MP_x_treat", "p"])
                interaction_hr = float(summ.loc["MP_x_treat", "exp(coef)"])
                main_treat_hr = float(summ.loc["treatment", "exp(coef)"])
                main_mp_hr = float(summ.loc["MP_z", "exp(coef)"])
            except Exception as e:
                C.log(f"  {sub_name}-{s}: cox err {e}")
                interaction_p = interaction_hr = main_treat_hr = main_mp_hr = np.nan

            # Median-split: ACT effect in high vs low MP
            med = df[s].median()
            df["MP_high"] = (df[s] >= med).astype(int)
            hr_act_high = hr_act_low = np.nan
            try:
                hi = df[df["MP_high"] == 1]
                if hi["treatment"].nunique() == 2:
                    cph_h = CoxPHFitter(penalizer=0.001)
                    cph_h.fit(hi[["treatment", "os_time", "os_event"]],
                                duration_col="os_time", event_col="os_event",
                                show_progress=False)
                    hr_act_high = float(cph_h.summary.loc["treatment", "exp(coef)"])
                lo = df[df["MP_high"] == 0]
                if lo["treatment"].nunique() == 2:
                    cph_l = CoxPHFitter(penalizer=0.001)
                    cph_l.fit(lo[["treatment", "os_time", "os_event"]],
                                duration_col="os_time", event_col="os_event",
                                show_progress=False)
                    hr_act_low = float(cph_l.summary.loc["treatment", "exp(coef)"])
            except Exception as e:
                C.log(f"  {sub_name}-{s}: split err {e}")

            rows.append({
                "subset": sub_name, "score": s, "n": len(df),
                "n_act": int((df["treatment"]==1).sum()),
                "n_obs": int((df["treatment"]==0).sum()),
                "main_mp_hr": main_mp_hr,
                "main_treat_hr": main_treat_hr,
                "interaction_hr": interaction_hr,
                "interaction_p": interaction_p,
                "hr_act_in_mp_high": hr_act_high,
                "hr_act_in_mp_low": hr_act_low,
            })

            for samp in df.index:
                box_rows.append({"subset": sub_name, "score": s,
                                  "Sample": samp,
                                  "value": float(df.loc[samp, s]),
                                  "treatment": "ACT" if df.loc[samp, "treatment"]==1 else "OBS",
                                  "os_time": float(df.loc[samp, "os_time"]),
                                  "os_event": int(df.loc[samp, "os_event"]),
                                  "MP_high": int(df.loc[samp, "MP_high"])})

    res = pd.DataFrame(rows)
    res.to_csv(C.OUT / f"{GSE}_response_comparison.csv", index=False)
    pd.DataFrame(box_rows).to_csv(C.OUT / f"{GSE}_boxplot_data.csv", index=False)
    C.log("\n=== Cox interaction MP × ACT (HR < 1 = ACT benefit; "
          "HR_high/HR_low ratio < 1 = MP-high benefits more) ===")
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
