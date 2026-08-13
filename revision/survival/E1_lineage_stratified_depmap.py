"""E1 — Lineage-stratified DepMap dependency (addresses non-LUAD specificity).

Main text compared 53 LUAD vs 1,097 mixed non-LUAD lines, ΔChronos = −0.077.
That difference is very small and does not support a claim of selectivity, and
pooling cancer types can in any case mask lineage specificity.

Three more informative alternatives:
  A lineage stratification: group by OncotreeLineage; rank mean LUAD Chronos among all lineages.
    Not specific and "most dependent among all lineages" are different claims; the latter still has value.
  B expression–dependency interaction: whether high-SEC61G LUAD lines are more dependent (mentioned in text, not in Methods).
  C comparison to established LUAD targets: place SEC61G lineage rank next to EGFR, KRAS, etc.,
    as an interpretable reference frame.
"""
import numpy as np, pandas as pd, os, warnings
from scipy import stats
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/depmap_lineage"
os.makedirs(OUT, exist_ok=True)
D = "${DATA_ROOT}/depmap/24Q2/"
GENES = ["SEC61G", "SRSF9", "ANGPTL4", "EGFR", "KRAS", "SEC61A1", "SEC61B"]

mod = pd.read_csv(D + "Model.csv", low_memory=False)
cr = pd.read_csv(D + "CRISPRGeneEffect.csv", index_col=0)
cr.columns = [c.split(" (")[0] for c in cr.columns]
cr = cr.loc[:, ~cr.columns.duplicated()]
m = mod.set_index("ModelID")[["OncotreeLineage", "OncotreeCode", "OncotreePrimaryDisease"]]
df = cr.join(m, how="inner")
print(f"cell lines {len(df)}, lineages {df.OncotreeLineage.nunique()}\n", flush=True)

print("=== A Lineage stratification: mean Chronos per lineage (more negative = more dependent) ===", flush=True)
for g in ["SEC61G", "EGFR", "KRAS"]:
    s = df.groupby("OncotreeLineage")[g].agg(n="size", mean="mean").query("n >= 15")
    s = s.sort_values("mean")
    luad_lines = df[df.OncotreeCode == "LUAD"]
    luad_mean = luad_lines[g].mean()
    rank_all = int((s["mean"] < luad_mean).sum()) + 1
    print(f"\n  {g}: LUAD ({len(luad_lines)} lines) mean {luad_mean:+.3f}")
    print(f"rank {rank_all} among {len(s)} lineages (n>=15); lower rank = more dependent")
    print(f"top 5 most dependent lineages:" +
          ", ".join(f"{i}({v['mean']:+.2f})" for i, v in s.head(5).iterrows()))
    print(f"Lung lineage overall {df[df.OncotreeLineage=='Lung'][g].mean():+.3f}"
          f"({int((df.OncotreeLineage=='Lung').sum())} lines)")
    s.to_csv(f"{OUT}/lineage_{g}.csv")

print("\n=== A2 LUAD vs each other lineage one-by-one (SEC61G) ===", flush=True)
luad_v = df.loc[df.OncotreeCode == "LUAD", "SEC61G"].dropna()
rows = []
for ln, s in df.groupby("OncotreeLineage"):
    v = s.loc[s.OncotreeCode != "LUAD", "SEC61G"].dropna()
    if len(v) < 15: continue
    u, p = stats.mannwhitneyu(luad_v, v, alternative="less")
    rows.append(dict(lineage=ln, n=len(v), mean=v.mean(),
                     delta=luad_v.mean() - v.mean(), p=p))
R = pd.DataFrame(rows).sort_values("delta")
R.to_csv(f"{OUT}/SEC61G_luad_vs_lineages.csv", index=False)
sig = R[(R.p < 0.05) & (R.delta < 0)]
print(f"lineages where LUAD is significantly more dependent: {len(sig)}/{len(R)}", flush=True)
print(R.head(6).round(4).to_string(index=False), flush=True)

print("\n=== B Expression–dependency interaction (within LUAD lines) ===", flush=True)
ex = pd.read_csv(D + "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", index_col=0,
                 usecols=lambda c: c == "Unnamed: 0" or c.split(" (")[0] in GENES)
ex.columns = [c.split(" (")[0] for c in ex.columns]
j = df.join(ex, how="inner", rsuffix="_expr")
lu = j[j.OncotreeCode == "LUAD"]
print(f"LUAD lines with expression + dependency data: {len(lu)}", flush=True)
for g in ["SEC61G", "SRSF9", "ANGPTL4"]:
    if g + "_expr" not in lu.columns and g not in ex.columns: continue
    e = lu[g + "_expr"] if g + "_expr" in lu.columns else lu[g]
    c = lu[g]
    ok = e.notna() & c.notna()
    if ok.sum() < 10: continue
    r, p = stats.spearmanr(e[ok], c[ok])
    ra, pa = stats.spearmanr(j[g + "_expr"] if g + "_expr" in j.columns else j[g],
                             j[g], nan_policy="omit")
    print(f"{g:<9} within LUAD, expression vs Chronos: ρ={r:+.3f} (p={p:.3g}, n={int(ok.sum())})"
          f"pan-cancer ρ={ra:+.3f} (p={pa:.2g})", flush=True)

print("\n=== C Reference: dependency strength and lineage rank per gene in LUAD ===", flush=True)
print(f"{'gene':<9}{'LUAD mean':>11}{'non-LUAD':>10}{'Δ':>9}{'one-sided p':>10}{'lineage rank':>10}", flush=True)
for g in GENES:
    if g not in df.columns: continue
    a = df.loc[df.OncotreeCode == "LUAD", g].dropna()
    b = df.loc[df.OncotreeCode != "LUAD", g].dropna()
    s = df.groupby("OncotreeLineage")[g].agg(n="size", mean="mean").query("n>=15").sort_values("mean")
    rk = int((s["mean"] < a.mean()).sum()) + 1
    print(f"  {g:<9}{a.mean():>+11.3f}{b.mean():>+10.3f}{a.mean()-b.mean():>+9.3f}"
          f"{stats.mannwhitneyu(a,b,alternative='less').pvalue:>10.2g}{f'{rk}/{len(s)}':>10}", flush=True)
print(f"\nResults written to {OUT}/", flush=True)
