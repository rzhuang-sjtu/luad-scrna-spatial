"""Step 24e: GSE42127 (UT Lung SPORE, n=176; 133 ADC) — Cox MP × ACT interaction.

ACT (mainly Carboplatin/Taxanes) vs OBS arms.
Same statistical framework as 24d.
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

GSE = "gse42127"


def main():
    C.log(f"=== {GSE.upper()} UT Lung SPORE ACT vs OBS ===")
    import GEOparse
    gse = GEOparse.get_GEO(geo="GSE42127", destdir=str(C.EXT), silent=True,
                             annotate_gpl=True, include_data=True)
    C.log(f"n_samples: {len(gse.gsms)}")

    expr = pd.DataFrame({g: s.table.set_index("ID_REF")["VALUE"]
                         for g, s in gse.gsms.items()})
    expr = expr.astype("float32")
    C.log(f"probe matrix: {expr.shape}")

    gpl = list(gse.gpls.values())[0]
    p2g = gpl.table[["ID", "Symbol"]].copy()
    p2g.columns = ["ID_REF", "gene"]
    p2g["gene"] = p2g["gene"].astype(str).str.strip()
    p2g = p2g[p2g["gene"].notna() & (p2g["gene"] != "") & (p2g["gene"] != "nan")]
    p2g = p2g.set_index("ID_REF")

    expr2 = expr.join(p2g, how="inner")
    expr_gene = expr2.groupby("gene")[expr.columns.tolist()].max()
    C.log(f"gene matrix: {expr_gene.shape}")

    gene_sets, ov = C.build_gene_sets(expr_gene.index)
    for k, n in ov.items():
        C.log(f"  {k}: {n} present")
    C.log("running ssGSEA")
    scores = C.run_ssgsea(expr_gene, gene_sets)
    C.log(f"  scores: {scores.shape}")

    # clinical
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
    C.log(f"adjuvant chemo: {clin['had_adjuvant_chemo'].value_counts().to_dict()}")
    C.log(f"histology: {clin['histology'].value_counts().to_dict()}")

    clin["os_time"] = pd.to_numeric(clin["overall_survival_months"], errors="coerce")
    clin["os_event"] = clin["survival_status"].astype(str).str.upper().eq("D").astype(int)
    clin["treatment"] = clin["had_adjuvant_chemo"].map(
        {"TRUE": 1, "FALSE": 0, "True": 1, "False": 0,
         True: 1, False: 0, "true": 1, "false": 0}).astype(float)

    score_df = scores.join(clin[["treatment", "os_time", "os_event",
                                    "histology", "final.pat.stage",
                                    "gender", "age_at_surgery"]],
                              how="inner")
    score_df = score_df.dropna(subset=["treatment", "os_time", "os_event"])
    score_df.index.name = "Sample"
    score_df.to_csv(C.OUT / f"{GSE}_mp_scores.csv")
    C.log(f"after QC: {score_df.shape}")

    from lifelines import CoxPHFitter
    SCORES = ["MP1", "MP2", "MP3", "MP4", "EMT_Hallmark",
               "Neutrophil_core", "NETs_composite"]
    SUBSETS = {
        "all_nsclc": score_df.index.tolist(),
        "ADC_only": score_df[score_df["histology"] == "Adenocarcinoma"].index.tolist(),
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
    C.log("\n=== Cox interaction MP × ACT ===")
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
