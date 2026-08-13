"""Compare 100-cell baseline (perturb/) vs 500-cell rerun (perturb_500/) Geneformer KO outputs.

Same filter as 04_venn.py: Sig==1 & Shift_to_goal_end>0 & FDR<0.05 & N_Detections>=5
Reports:
  - per-transition top-50/100/200 overlap (Jaccard, kept/dropped/new genes)
  - rank stability of original 7-gene candidate pool in the 500-cell run
  - new 500-cell ≥2-hit pool vs original
"""
from pathlib import Path
import pandas as pd

ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer")
P100, P500 = ROOT / "perturb", ROOT / "perturb_500"
OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn_500_diff")
OUT.mkdir(parents=True, exist_ok=True)

TRANSITIONS = ["macro_spp1_to_c1qc", "mal_mp3_to_mp1", "neu_osm_priming_to_low"]
N_MIN, FDR = 5, 0.05
TOP_NS = [50, 100, 200]

ORIGINAL_POOL = ["CNBP", "HSP90AB1", "NDUFB2", "PARK7", "SLC25A5", "SRSF9", "TMEM59"]


def load_filtered(base, t):
    df = pd.read_csv(base / t / f"{t}_stats.csv", index_col=0)
    df = df[(df["Sig"] == 1) & (df["Shift_to_goal_end"] > 0)
            & (df["Goal_end_FDR"] < FDR) & (df["N_Detections"] >= N_MIN)].copy()
    return df.sort_values("Shift_to_goal_end", ascending=False).reset_index(drop=True)


def gset(d, n):
    return set(d.head(n)["Gene_name"])


def jaccard(a, b):
    if not (a or b): return float("nan")
    return len(a & b) / len(a | b)


print("=" * 78)
print(" 100 vs 500 cell comparison — Geneformer KO Fig 8 candidates")
print("=" * 78)

# Per-transition pool sizes
print("\n[A] Pool sizes after filter (Sig=1, FDR<0.05, ΔS>0, n_det≥5)")
print(f"  {'transition':<28}  {'100-cell':>10}  {'500-cell':>10}")
size_rows = []
filt100, filt500 = {}, {}
for t in TRANSITIONS:
    filt100[t] = load_filtered(P100, t)
    filt500[t] = load_filtered(P500, t)
    print(f"  {t:<28}  {len(filt100[t]):>10}  {len(filt500[t]):>10}")
    size_rows.append({"transition": t, "n_100": len(filt100[t]), "n_500": len(filt500[t])})
pd.DataFrame(size_rows).to_csv(OUT / "pool_sizes.csv", index=False)

# Per-transition top-N overlap
print("\n[B] Top-N overlap per transition (Jaccard | kept | dropped | new)")
overlap_rows = []
for t in TRANSITIONS:
    print(f"\n  {t}")
    for n in TOP_NS:
        a = gset(filt100[t], n); b = gset(filt500[t], n)
        kept = a & b; dropped = a - b; new = b - a
        print(f"    top-{n:>3}: J={jaccard(a, b):.3f}  kept={len(kept)}  dropped={len(dropped)}  new={len(new)}")
        overlap_rows.append({
            "transition": t, "top_N": n,
            "jaccard": jaccard(a, b),
            "n_kept": len(kept), "n_dropped": len(dropped), "n_new": len(new),
            "dropped_genes": ";".join(sorted(dropped)),
            "new_genes": ";".join(sorted(new)),
        })
pd.DataFrame(overlap_rows).to_csv(OUT / "topN_overlap_by_transition.csv", index=False)

# Original 7-gene pool fate in 500-cell at top-200 per transition
print("\n[C] Original 7-gene ≥2-hit pool fate in 500-cell run (top-200 per transition)")
print(f"  {'gene':<10}  {'macro':>10}  {'mal':>10}  {'neu':>10}  {'n_hit_500':>10}  status")
pool_fate_rows = []
top200_500 = {t: gset(filt500[t], 200) for t in TRANSITIONS}
top200_100 = {t: gset(filt100[t], 200) for t in TRANSITIONS}
for g in ORIGINAL_POOL:
    ranks_500, ranks_100 = {}, {}
    for t in TRANSITIONS:
        d100, d500 = filt100[t], filt500[t]
        h100 = d100.index[d100["Gene_name"] == g]
        h500 = d500.index[d500["Gene_name"] == g]
        ranks_100[t] = int(h100[0]) + 1 if len(h100) else None
        ranks_500[t] = int(h500[0]) + 1 if len(h500) else None
    n_hit_500 = sum(1 for t in TRANSITIONS if g in top200_500[t])
    n_hit_100 = sum(1 for t in TRANSITIONS if g in top200_100[t])
    status = "PRESERVED ≥2-hit" if n_hit_500 >= 2 else ("DROPPED" if n_hit_100 >= 2 else "n/a")
    def fmt(r): return f"{r}" if r is not None else "—"
    print(f"  {g:<10}  {fmt(ranks_500['macro_spp1_to_c1qc']):>10}  "
          f"{fmt(ranks_500['mal_mp3_to_mp1']):>10}  "
          f"{fmt(ranks_500['neu_osm_priming_to_low']):>10}  "
          f"{n_hit_500:>10}  {status}")
    pool_fate_rows.append({
        "gene": g,
        "rank_macro_100": ranks_100["macro_spp1_to_c1qc"], "rank_macro_500": ranks_500["macro_spp1_to_c1qc"],
        "rank_mal_100":   ranks_100["mal_mp3_to_mp1"],     "rank_mal_500":   ranks_500["mal_mp3_to_mp1"],
        "rank_neu_100":   ranks_100["neu_osm_priming_to_low"], "rank_neu_500": ranks_500["neu_osm_priming_to_low"],
        "n_hit_in_top200_100": n_hit_100,
        "n_hit_in_top200_500": n_hit_500,
        "status": status,
    })
pd.DataFrame(pool_fate_rows).to_csv(OUT / "original_pool_fate.csv", index=False)

# New 500-cell ≥2-hit pool at top-200
print("\n[D] 500-cell ≥2-hit pool at top-200 (compare to original 7)")
a, b, c = top200_500["macro_spp1_to_c1qc"], top200_500["mal_mp3_to_mp1"], top200_500["neu_osm_priming_to_low"]
pool_500 = (a & b) | (a & c) | (b & c)
pool_500_3way = a & b & c
pool_100 = set(ORIGINAL_POOL)
print(f"  500-cell ≥2-hit pool size: {len(pool_500)} (vs 100-cell: {len(pool_100)})")
print(f"    preserved (in both):    {sorted(pool_100 & pool_500)}")
print(f"    lost (only in 100):     {sorted(pool_100 - pool_500)}")
print(f"    gained (only in 500):   {sorted(pool_500 - pool_100)}")
print(f"  500-cell 3-way intersection: {sorted(pool_500_3way) or '∅'}")

pd.DataFrame([{"gene": g, "in_100_pool": g in pool_100, "in_500_pool": g in pool_500,
               "n_hits_500": sum(1 for s in (a, b, c) if g in s)}
              for g in sorted(pool_100 | pool_500)]
             ).to_csv(OUT / "pool_diff_summary.csv", index=False)

# write 500-cell pool CSV with full per-transition stats (mirror 8A_candidate_pool.csv format)
pool_500_rows = []
for g in sorted(pool_500):
    row = {"Gene_name": g}; hits = []
    for t in TRANSITIONS:
        d = filt500[t]
        sub = d[d["Gene_name"] == g]
        if len(sub):
            short = t.split("_")[0]
            row[f"{short}_rank"]  = int(d.index[d["Gene_name"] == g][0]) + 1
            row[f"{short}_shift"] = sub["Shift_to_goal_end"].iloc[0]
            row[f"{short}_FDR"]   = sub["Goal_end_FDR"].iloc[0]
            hits.append(t)
            if "Ensembl_ID" not in row: row["Ensembl_ID"] = sub["Ensembl_ID"].iloc[0]
    row["n_transitions_hit"] = len(hits)
    row["transitions_hit"] = ";".join(hits)
    pool_500_rows.append(row)
pd.DataFrame(pool_500_rows).sort_values(["n_transitions_hit", "Gene_name"], ascending=[False, True]
    ).to_csv(OUT / "8A_500_candidate_pool.csv", index=False)

print(f"\n[E] CSV outputs in: {OUT}")
print("    pool_sizes.csv")
print("    topN_overlap_by_transition.csv")
print("    original_pool_fate.csv")
print("    pool_diff_summary.csv")
print("    8A_500_candidate_pool.csv")
print("\nDONE.")
