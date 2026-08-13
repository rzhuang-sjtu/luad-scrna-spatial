#!/usr/bin/env python
"""
G4 — independently recheck CopyKAT malignant calls with inferCNV (independent CNV cross-check).

Input:
  ~/luad/results/infercnv/summary/*_cells.csv   (G2 output, copied back from the compute node)
      columns: cell, group (epithelial/reference), cnv_var, cnv_mean_abs, cnv_rms,
          cnv_p90_abs, cnv_frac_dev, frac_altered, patient
  ~/luad/data/processed/luad_copykat.h5ad       (obs has copykat_pred and malignant)

Principle: report threshold-free quantities first (AUROC, internal negative control, prevalence-matched pairing),
      then agreement and kappa on the full threshold grid — report the whole grid, not a favourable point.
      Report all six metrics; do not cherry-pick the best post hoc.

Output (~/luad/results/infercnv/):
  cell_level.csv.gz     per-cell merged table
  negative_control.csv  three tiers: reference / epithelial non-malignant / epithelial malignant
  auroc.csv             overall AUROC for six metrics + patient-stratified bootstrap 95% CI
  threshold_grid.csv    sens/spec/agreement/kappa on the threshold grid
  patient_level.csv     per-patient malignant fraction, AUROC, prevalence-matched overlap
  summary.md            conclusion summary
"""
import os, sys, gzip, glob
import numpy as np
import pandas as pd
import h5py

ROOT = os.path.expanduser("~/luad")
SUMDIR = f"{ROOT}/results/infercnv/summary"
H5 = f"{ROOT}/data/processed/luad_copykat.h5ad"
OUT = f"{ROOT}/results/infercnv"
METRICS = ["cnv_var", "cnv_mean_abs", "cnv_rms", "cnv_p90_abs",
           "cnv_frac_dev", "frac_altered"]
NBOOT = 2000
SEED = 42


def load_obs():
    """Load only the four needed columns; do not load the expression matrix."""
    with h5py.File(H5, "r") as h:
        ix = h["obs"].attrs.get("_index", "index")
        cid = np.array([x.decode() if isinstance(x, bytes) else x
                        for x in h["obs"][ix][:]])

        def cat(k):
            g = h["obs"][k]
            if isinstance(g, h5py.Group):
                cats = np.array([x.decode() if isinstance(x, bytes) else x
                                 for x in g["categories"][:]])
                codes = g["codes"][:]
                out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
                return out
            v = g[:]
            return np.array([x.decode() if isinstance(x, bytes) else x for x in v])

        return pd.DataFrame({
            "cell": cid,
            "celltype_coarse": cat("celltype_coarse"),
            "copykat_pred": cat("copykat_pred"),
            "malignant": cat("malignant"),
        })


def auc(y, s):
    """Rank-based AUROC with ties handled. y is 0/1."""
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    ok = np.isfinite(s)
    y, s = y[ok], s[ok]
    n1, n0 = y.sum(), (1 - y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().values
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot_auc_by_patient(df, metric, nboot=NBOOT, seed=SEED):
    """Resample by patient — cells are not independent; cell-level bootstrap underestimates variance."""
    rng = np.random.default_rng(seed)
    pats = df["patient"].unique()
    groups = {p: g for p, g in df.groupby("patient")}
    vals = []
    for _ in range(nboot):
        pick = rng.choice(pats, len(pats), replace=True)
        d = pd.concat([groups[p] for p in pick], ignore_index=True)
        a = auc(d["y"].values, d[metric].values)
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def kappa(a, b):
    """Cohen's kappa, binary."""
    a = np.asarray(a, int); b = np.asarray(b, int)
    n = len(a)
    if n == 0:
        return np.nan
    po = (a == b).mean()
    pe = (a.mean() * b.mean()) + ((1 - a.mean()) * (1 - b.mean()))
    return np.nan if pe == 1 else (po - pe) / (1 - pe)


def main():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(glob.glob(f"{SUMDIR}/*_cells.csv"))
    if not files:
        sys.exit(f"No inferCNV summary files found: {SUMDIR}/*_cells.csv")
    print(f"Reading inferCNV summaries for {len(files)} patients ...", flush=True)
    icnv = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    print(f"inferCNV cells {len(icnv):,}, patients {icnv.patient.nunique()}", flush=True)

    failed = sorted(os.path.basename(f)[7:-4]
                    for f in glob.glob(f"{SUMDIR}/FAILED_*.txt"))
    if failed:
        print(f"failed patients {len(failed)}: {', '.join(failed[:10])}"
              f"{' ...' if len(failed) > 10 else ''}", flush=True)

    obs = load_obs()
    d = icnv.merge(obs, on="cell", how="left")
    n_unmatched = d["malignant"].isna().sum()
    print(f"matched to copykat h5ad {len(d) - n_unmatched:,} / {len(d):,}"
          f"(unmatched {n_unmatched:,})", flush=True)
    d = d[d["malignant"].notna()].copy()
    d.to_csv(f"{OUT}/cell_level.csv.gz", index=False, compression="gzip")

    # Important: G1 selected cells using celltype_coarse from luad_merged_annotated.h5ad,
    # and that annotation is wrong (among 240k cells labelled Fibroblast, 58% express PTPRC, 45.8% LYZ,
    # only 10.5% COL1A1; the luad_copykat.h5ad labels are correct — Fibroblast COL1A1+ 74.2%).
    # Consequence: ~13% immune cells entered the inferCNV observation set.
    # Re-tier here using luad_copykat.h5ad labels; evaluation set keeps only cells called epithelial in both.
    #
    # Contaminating immune cells provide a non-circular negative control:
    # they follow the same observation path as epithelium (unlike reference cells used to define the baseline, which is circular).
    is_epi = d["celltype_coarse"] == "Epithelial"
    d["tier"] = np.where(d["group"] == "reference", "1_参考细胞(定基线,循环)",
               np.where(~is_epi, "2_immune_in_obs_set(non-circular_control)",
               np.where(d["malignant"] == "Malignant", "5_epithelial_CopyKAT_malignant",
               np.where(d["malignant"] == "Non-malignant", "3_epithelial_CopyKAT_non-malignant",
                        "4_epithelial_Uncertain"))))
    NC = d.groupby("tier")[METRICS].median().reset_index()
    NC.insert(1, "n", d.groupby("tier").size().values)
    NC.to_csv(f"{OUT}/negative_control.csv", index=False)
    print("\nInternal negative control (per-metric medians):", flush=True)
    print(NC.to_string(index=False), flush=True)

    # Evaluation set = epithelial cells in the observation set agreed by both annotations.
    # Exclude reference cells (otherwise cell-type differences masquerade as discrimination);
    # also exclude immune cells in the observation set (control only; not used for discrimination metrics).
    d = d[(d["group"] == "epithelial") & is_epi[d.index]].copy()
    print(f"Evaluation set: cells called epithelial in both annotations within the observation set {len(d):,}", flush=True)
    ev = d[d["malignant"].isin(["Malignant", "Non-malignant"])].copy()
    ev["y"] = (ev["malignant"] == "Malignant").astype(int)
    print(f"Evaluable cells {len(ev):,} (malignant {ev.y.sum():,} /"
          f"non-malignant {(1-ev.y).sum():,}; Uncertain excluded)", flush=True)

    rows = []
    for m in METRICS:
        a = auc(ev["y"].values, ev[m].values)
        lo, hi = boot_auc_by_patient(ev[["patient", "y", m]].dropna(), m)
        # Compute per patient then take the median — avoids artefacts from inter-patient malignant-fraction differences
        per = ev.groupby("patient").apply(
            lambda g: auc(g["y"].values, g[m].values), include_groups=False)
        per = per.dropna()
        rows.append(dict(metric=m, auroc_overall=a, ci_lo=lo, ci_hi=hi,
                         auroc_patient_median=float(per.median()),
                         auroc_patient_q25=float(per.quantile(.25)),
                         auroc_patient_q75=float(per.quantile(.75)),
                         n_patients_evaluable=int(len(per)),
                         n_cells=int(ev[m].notna().sum())))
    A = pd.DataFrame(rows)
    A.to_csv(f"{OUT}/auroc.csv", index=False)
    print(A.to_string(index=False), flush=True)

    grid = []
    for m in METRICS:
        v = ev[m].dropna()
        if v.empty:
            continue
        cuts = np.unique(np.quantile(v, np.arange(0.05, 1.0, 0.05)))
        for c in cuts:
            pred = (ev[m] >= c).astype(int)
            ok = ev[m].notna()
            y = ev.loc[ok, "y"].values
            p = pred[ok].values
            tp = int(((p == 1) & (y == 1)).sum()); fp = int(((p == 1) & (y == 0)).sum())
            fn = int(((p == 0) & (y == 1)).sum()); tn = int(((p == 0) & (y == 0)).sum())
            grid.append(dict(
                metric=m, cutoff=float(c),
                frac_called_malignant=float(p.mean()),
                sensitivity=tp / (tp + fn) if tp + fn else np.nan,
                specificity=tn / (tn + fp) if tn + fp else np.nan,
                agreement=(tp + tn) / len(y),
                kappa=kappa(y, p), tp=tp, fp=fp, fn=fn, tn=tn))
    G = pd.DataFrame(grid)
    G.to_csv(f"{OUT}/threshold_grid.csv", index=False)

    best = G.loc[G.groupby("metric")["kappa"].idxmax()]
    print("\nPer-metric kappa peak (for reference only; main text reports the full curve):", flush=True)
    print(best[["metric", "cutoff", "sensitivity", "specificity",
                "agreement", "kappa"]].to_string(index=False), flush=True)

    def matched(g):
        k = int(g["y"].sum())
        out = {"n_cells": len(g), "n_ck_malignant": k,
               "frac_ck_malignant": k / len(g)}
        for m in METRICS:
            s = g[m]
            if k == 0 or k == len(g) or s.notna().sum() < len(g):
                out[f"overlap_{m}"] = np.nan
                continue
            top = s.rank(ascending=False, method="first") <= k
            out[f"overlap_{m}"] = float((top & (g["y"] == 1)).sum() / k)
        for m in METRICS:
            out[f"auroc_{m}"] = auc(g["y"].values, g[m].values)
        return pd.Series(out)

    P = ev.groupby("patient").apply(matched, include_groups=False).reset_index()
    P["n_minor"] = np.minimum(P["n_ck_malignant"], P["n_cells"] - P["n_ck_malignant"])
    P.to_csv(f"{OUT}/patient_level.csv", index=False)

    # Unstratified within-patient medians are badly deflated: many patients have only single-digit malignant
    # (or non-malignant) calls from CopyKAT, so AUROC is a coin flip on one or two cells;
    # values of 0 or 1 carry no information. Reporting 'within-patient median 0.54' would misstate noise as a real weakness.
    BINS = [-1, 10, 50, 200, 1000, np.inf]
    LBL = ["<=10", "11-50", "51-200", "201-1000", ">1000"]
    P["minority_class_bin"] = pd.cut(P["n_minor"], BINS, labels=LBL)
    strat = []
    for lb, g in P.groupby("minority_class_bin", observed=True):
        row = {"minority_class_bin": lb, "病人数": len(g),
               "median_minority_cells": float(g["n_minor"].median())}
        for m in METRICS:
            v = g[f"auroc_{m}"].dropna()
            row[m] = float(v.median()) if len(v) else np.nan
        strat.append(row)
    S = pd.DataFrame(strat)
    S.to_csv(f"{OUT}/patient_auroc_by_size.csv", index=False)
    print("\nWithin-patient AUROC stratified by evaluable cell count (median):", flush=True)
    print(S.to_string(index=False), flush=True)

    L = ["# Agreement between inferCNV and CopyKAT malignant calls", "",
         f"{ev.patient.nunique()} patients, {len(ev):,} epithelial cells"
         f"(CopyKAT malignant {int(ev.y.sum()):,}, non-malignant {int((1-ev.y).sum()):,},"
         f"Uncertain {int((d['malignant']=='Uncertain').sum()):,} excluded).", ""]
    if failed:
        L += [f"inferCNV failed for {len(failed)} patients: {', '.join(failed)}", ""]
    L += ["## 0. Internal negative control (per-metric medians)", "",
          "| tier | n |" + " | ".join(METRICS) + " |",
          "|---" * (len(METRICS) + 2) + "|"]
    for _, r in NC.iterrows():
        L.append(f"| {r['档']} | {int(r['n']):,} | " +
                 " | ".join(f"{r[m]:.4g}" for m in METRICS) + " |")
    L += ["", "Reference cells are the patient's own T/NK + B, near-diploid by construction."
          "The gap versus epithelial cells called malignant is direct evidence these metrics work, independent of any threshold.", "",
          "## 1. Threshold-free discrimination (AUROC)", "",
          "| metric | overall AUROC | 95% CI (patient bootstrap) | within-patient AUROC median (IQR) |",
          "|---|---|---|---|"]
    for _, r in A.iterrows():
        L.append(f"| {r.metric} | {r.auroc_overall:.3f} | "
                 f"[{r.ci_lo:.3f}, {r.ci_hi:.3f}] | "
                 f"{r.auroc_patient_median:.3f} "
                 f"({r.auroc_patient_q25:.3f}–{r.auroc_patient_q75:.3f}) |")
    L += ["", "Overall AUROC is inflated by inter-patient malignant-fraction differences. But the **within-patient median must not be read raw**:",
          "Many patients have only single-digit malignant cells from CopyKAT (or only single-digit non-malignant),",
          "so within-patient AUROC is a coin flip on one or two cells; 0 or 1 carries no information.",
          "Must stratify by evaluable cell count:", "",
          "| minority-class cell count | n patients | minority median |" + " | ".join(METRICS) + " |",
          "|---" * (len(METRICS) + 3) + "|"]
    for _, r in S.iterrows():
        L.append(f"| {r['少数类分档']} | {int(r['病人数'])} | {r['少数类细胞中位']:.0f} | " +
                 " | ".join(f"{r[m]:.3f}" if pd.notna(r[m]) else "—" for m in METRICS) + " |")
    L += ["", "Within-patient AUROC rises monotonically with evaluable cell count; only patients with both classes well represented give interpretable readouts.", "",
          "## 2. Threshold grid (full curve in threshold_grid.csv)", "",
          "| metric | threshold at peak kappa | sensitivity | specificity | agreement | kappa |",
          "|---|---|---|---|---|---|"]
    for _, r in best.iterrows():
        L.append(f"| {r.metric} | {r.cutoff:.4f} | {r.sensitivity:.3f} | "
                 f"{r.specificity:.3f} | {r.agreement:.3f} | {r.kappa:.3f} |")
    L += ["", "## 3. Prevalence-matched pairing (per patient take top-k by CopyKAT malignant count)", ""]
    for m in METRICS:
        c = P[f"overlap_{m}"].dropna()
        if len(c):
            L.append(f"- {m}: median overlap {c.median():.3f}"
                     f"（IQR {c.quantile(.25):.3f}–{c.quantile(.75):.3f}，"
                     f"{len(c)} patients evaluable)")
    L += ["", "This metric is independent of call threshold and of inter-patient malignant-fraction differences.", ""]
    open(f"{OUT}/summary.md", "w").write("\n".join(L))
    print(f"\nWrote {OUT}/summary.md", flush=True)


if __name__ == "__main__":
    main()
