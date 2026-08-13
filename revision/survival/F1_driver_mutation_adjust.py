"""F1 — MP3/MP4 survival analysis adjusted for EGFR/KRAS driver mutations.

Background: the main-text limitations said that adjustment for driver mutations
(EGFR/KRAS) was limited by missing annotation and would need larger cohorts.
The MC3 non-silent mutation matrix (mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz) is
available, so that caveat no longer holds and can be answered with data.

Three tasks:
  1. Frequencies of EGFR / KRAS / TP53 / STK11 / KEAP1 mutations in TCGA-LUAD, checked against literature
  2. Whether each MP score differs between mutant and wild-type (potential confounding)
  3. Add mutation status to the Cox model and check stability of MP3 / MP4 / composite-score HRs
"""
import numpy as np, pandas as pd, os, warnings
from lifelines import CoxPHFitter
from scipy import stats
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/driver_mutation"
os.makedirs(OUT, exist_ok=True)
DRIVERS = ["EGFR", "KRAS", "TP53", "STK11", "KEAP1"]

mc3 = pd.read_csv("${WORK_ROOT}/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz",
                  sep="\t", index_col=0)
print(f"MC3 matrix {mc3.shape} (genes x samples)", flush=True)
mc3.columns = [c[:15] for c in mc3.columns]
mut = mc3.loc[[g for g in DRIVERS if g in mc3.index]].T
mut.columns = [f"mut_{c}" for c in mut.columns]

S = pd.read_csv("${WORK_ROOT}/luad_figures/fig3/tcga_luad_mp_ssgsea.csv.gz", index_col=0)
S.index = [i[:15] for i in S.index]
clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin = clin[clin.sample_type == "Primary Tumor"].copy()
clin["s15"] = clin.sample_barcode.str[:15]
clin = clin.drop_duplicates("s15").set_index("s15")

d = clin.join(S, how="inner").join(mut, how="inner")
print(f"three-way intersecting samples {len(d)}", flush=True)
d["event"] = (d.vital_status.str.strip().str.lower() == "dead").astype(int)
d["time"] = np.where(d.event == 1, d.days_to_death, d.days_to_last_follow_up)
d = d[d.time.notna() & (d.time > 0)]
d["age"] = pd.to_numeric(d.age_at_diagnosis, errors="coerce") / 365.25
d["stage_num"] = d.ajcc_stage.astype(str).str.extract(r"(Stage [IV]+)", expand=False).map(
    {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
d["male"] = (d.gender.str.lower() == "male").astype(int)
for c in ["MP1", "MP2", "MP3", "MP4"]:
    d["z" + c] = (d[c] - d[c].mean()) / d[c].std()
d["dual"] = d.zMP3 - d.zMP4

print("\n=== 1 Driver-mutation frequencies (vs literature) ===", flush=True)
lit = {"EGFR": "~14% (Western LUAD)", "KRAS": "约 33%", "TP53": "约 46%",
       "STK11": "~17%", "KEAP1": "约 17%"}
for g in DRIVERS:
    c = f"mut_{g}"
    if c not in d.columns: continue
    print(f"{g:<7} {d[c].mean()*100:5.1f}%  ({int(d[c].sum())}/{len(d)})   literature {lit[g]}", flush=True)

print("\n=== 2 MP score differences mutant vs wild-type (potential confounding) ===", flush=True)
rows = []
for g in DRIVERS:
    c = f"mut_{g}"
    if c not in d.columns: continue
    line = f"  {g:<7}"
    for mp in ["MP1", "MP2", "MP3", "MP4"]:
        a, b = d.loc[d[c] == 1, "z" + mp], d.loc[d[c] == 0, "z" + mp]
        u, p = stats.mannwhitneyu(a, b)
        line += f"  {mp}:Δ={a.mean()-b.mean():+.2f}({'*' if p<0.05 else ' '})"
        rows.append(dict(driver=g, mp=mp, delta=a.mean()-b.mean(), p=p))
    print(line, flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/mp_by_driver.csv", index=False)

print("\n=== 3 Cox models: HR stability for MPs after adding driver mutations ===", flush=True)
base = ["age", "stage_num", "male"]
mutc = [f"mut_{g}" for g in DRIVERS if f"mut_{g}" in d.columns]
sets = [("baseline (age/stage/sex)", base),
        ("+ EGFR/KRAS", base + [c for c in mutc if c.endswith(("EGFR", "KRAS"))]),
        ("+ all five driver genes", base + mutc)]
for target, lab in [(["zMP1", "zMP2", "zMP3", "zMP4"], "all four MPs jointly"), (["dual"], "复合分")]:
    print(f"\n  【{lab}】", flush=True)
    for name, cov in sets:
        cols = ["time", "event"] + cov + target
        sub = d[cols].dropna()
        try:
            c_ = CoxPHFitter().fit(sub, "time", "event")
            out = "  ".join(
                f"{t.replace('z','')}: HR={c_.summary.loc[t,'exp(coef)']:.3f}"
                f"({c_.summary.loc[t,'exp(coef) lower 95%']:.2f}-{c_.summary.loc[t,'exp(coef) upper 95%']:.2f})"
                f" p={c_.summary.loc[t,'p']:.3g}" for t in target)
            print(f"{name:<22} n={len(sub)} events={int(sub.event.sum())}  {out}", flush=True)
        except Exception as e:
            print(f"{name}: fit failed {type(e).__name__}", flush=True)
print(f"\nResults written to {OUT}/", flush=True)
