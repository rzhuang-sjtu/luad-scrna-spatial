"""Step 20c (Fig 3I/J): Drug sensitivity per MP via DepMap PRISM ridge regression.

Pipeline:
  1. Download DepMap PRISM Repurposing AUC + cell-line expression (Sanger TPM).
     Cached under ${DATA_ROOT}/depmap/.
  2. Train ridge regression per drug: cell-line gene expression → AUC.
  3. Predict on MP pseudo-bulk (mean log1p expression per MP).
  4. Per-MP rank drugs; identify selective drugs (high variance across MPs).

Outputs to ${WORK_ROOT}/luad_figures/fig3/:
  drug_sensitivity_per_mp.csv     drug × MP predicted AUC
  drug_selectivity_summary.csv    per-MP top selective drugs
  drug_sensitivity_summary.md
"""
from __future__ import annotations
import os, time, urllib.request
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

CACHE = Path("${DATA_ROOT}/depmap")
CACHE.mkdir(parents=True, exist_ok=True)
MAL = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
FIG3 = Path("${WORK_ROOT}/luad_figures/fig3")

# DepMap 24Q4 release URLs (stable). PRISM Repurposing Secondary AUC.
URLS = {
    # Cell-line expression (TPM log2(x+1))
    "expression": ("https://figshare.com/ndownloader/files/40448555",
                    CACHE/"OmicsExpressionProteinCodingGenesTPMLogp1.csv"),
    # PRISM Repurposing primary AUC (cell line × drug)
    "prism": ("https://figshare.com/ndownloader/files/40448561",
               CACHE/"PRISM_Repurposing_Primary_Screen.csv"),
}


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 1_000_000:  # >1MB to consider real
        log(f"  cache hit: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return True
    log(f"  downloading {url} → {dest.name}")
    try:
        urllib.request.urlretrieve(url, dest)
        sz = dest.stat().st_size
        log(f"    fetched {sz/1e6:.2f} MB")
        if sz < 1_000_000:
            log(f"    too small ({sz} bytes) — treating as fail")
            dest.unlink(missing_ok=True)
            return False
        return True
    except Exception as e:
        log(f"    FAIL: {e}")
        return False


def fallback_msigdb_cgp_approach(mal_h5ad, fig3_dir):
    """If DepMap download fails, use MSigDB C2:CGP perturbation gene sets
    via gseapy.ssgsea on per-MP signatures."""
    log("FALLBACK: MSigDB C2:CGP perturbation enrichment")
    # Download MSigDB C2:CGP gmt
    cgp_path = Path.home()/"luad/data/gmt/c2.cgp.v2024.1.Hs.symbols.gmt"
    if not cgp_path.exists():
        cgp_url = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
                    "2024.1.Hs/c2.cgp.v2024.1.Hs.symbols.gmt")
        try:
            urllib.request.urlretrieve(cgp_url, cgp_path)
            log(f"  fetched CGP: {cgp_path.stat().st_size/1e6:.1f} MB")
        except Exception as e:
            log(f"  CGP download FAIL: {e}")
            return False

    # Build per-MP pseudo-bulk
    a = sc.read_h5ad(mal_h5ad)
    if a.X.min() < 0:
        a.X = np.log1p(a.layers["counts"].astype(np.float32))
    pb_rows = []
    mps = ["MP1","MP2","MP3","MP4"]
    for mp in mps:
        mask = (a.obs["dominant_MP"].astype(str) == mp).values
        m = np.asarray(a[mask].X.mean(axis=0)).ravel()
        pb_rows.append(m)
    pb = pd.DataFrame(np.array(pb_rows).T, index=a.var_names, columns=mps)
    log(f"  pseudo-bulk for ssGSEA: {pb.shape}")

    import gseapy as gp
    ss = gp.ssgsea(data=pb, gene_sets=str(cgp_path), outdir=None,
                    sample_norm_method="rank", no_plot=True,
                    min_size=10, max_size=2000, permutation_num=0,
                    seed=0, threads=8)
    res = ss.res2d.pivot_table(index="Term", columns="Name", values="NES").astype(float)
    log(f"  CGP × MP NES: {res.shape}")
    res.to_csv(fig3_dir/"drug_sensitivity_per_mp.csv")

    # Selectivity = std across MPs
    res["std_across_mp"] = res[mps].std(axis=1)
    res["mean_NES"] = res[mps].mean(axis=1)
    sel = res.sort_values("std_across_mp", ascending=False)
    # Top 20 selective per MP (where that MP is the max)
    rows = []
    for mp in mps:
        mp_top = res[mps].copy()
        mp_top["selective_for"] = res[mps].idxmax(axis=1)
        mp_for = mp_top[mp_top["selective_for"] == mp].sort_values(mp, ascending=False).head(20)
        for term, vals in mp_for.iterrows():
            rows.append({"MP": mp, "perturbation": term, "NES": vals[mp],
                          "std_across_mp": float(res.loc[term, "std_across_mp"])})
    sel_df = pd.DataFrame(rows)
    sel_df.to_csv(fig3_dir/"drug_selectivity_summary.csv", index=False)
    log(f"  selectivity summary saved ({len(sel_df)})")

    with open(fig3_dir/"drug_sensitivity_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 20c — Drug/perturbation sensitivity per MP\n\n")
        f.write("**METHOD**: MSigDB C2:CGP perturbation gene-set enrichment (fallback)\n")
        f.write("DepMap PRISM download failed; used MSigDB chemical/genetic perturbation\n"
                "gene sets via ssGSEA on per-MP pseudo-bulks. NES > 0 means perturbation\n"
                "signature is up-regulated in that MP.\n\n")
        f.write(f"- N gene sets: {len(res)}\n- Output: `drug_sensitivity_per_mp.csv` (Term × MP NES)\n")
        f.write(f"- Per-MP top-20 selective: `drug_selectivity_summary.csv`\n")
    return True


def main():
    t0 = time.time()
    # Try DepMap PRISM first
    log("attempting DepMap PRISM download")
    ok_expr = fetch(*URLS["expression"])
    ok_prism = fetch(*URLS["prism"])

    if not (ok_expr and ok_prism):
        log("DepMap fetch incomplete, using MSigDB CGP fallback")
        fallback_msigdb_cgp_approach(MAL, FIG3)
        log(f"DONE (fallback) in {time.time()-t0:.1f}s")
        return

    # Load DepMap data
    log("loading DepMap expression + PRISM AUC")
    expr = pd.read_csv(URLS["expression"][1], index_col=0)
    log(f"  expression: {expr.shape}")
    prism = pd.read_csv(URLS["prism"][1], index_col=0)
    log(f"  PRISM: {prism.shape}")

    # Align cell lines
    common = expr.index.intersection(prism.index)
    log(f"  common cell lines: {len(common)}")
    expr = expr.loc[common]
    prism = prism.loc[common]

    # Strip parens from gene names: "CD8A (925)" -> "CD8A"
    expr.columns = expr.columns.str.replace(r"\s*\(\d+\)$", "", regex=True)

    # Build per-MP pseudo-bulk
    log("building per-MP pseudo-bulks")
    a = sc.read_h5ad(MAL)
    if a.X.min() < 0:
        a.X = np.log1p(a.layers["counts"].astype(np.float32))
    pb_rows = []; mps = ["MP1","MP2","MP3","MP4"]
    for mp in mps:
        mask = (a.obs["dominant_MP"].astype(str) == mp).values
        m = np.asarray(a[mask].X.mean(axis=0)).ravel()
        pb_rows.append(m)
    pb = pd.DataFrame(np.array(pb_rows), index=mps, columns=a.var_names)
    log(f"  pseudo-bulk: {pb.shape}")

    # Match gene universe
    common_genes = sorted(set(expr.columns) & set(pb.columns))
    log(f"  common genes: {len(common_genes)}")
    X_train = expr[common_genes].values
    X_pred = pb[common_genes].values

    log("fitting ridge regression per drug")
    from sklearn.linear_model import Ridge
    drugs = prism.columns.tolist()
    pred = np.full((len(mps), len(drugs)), np.nan)
    n_ok = 0
    for j, drug in enumerate(drugs):
        y = prism[drug].values
        ok = np.isfinite(y) & np.all(np.isfinite(X_train), axis=1)
        if ok.sum() < 50: continue
        try:
            r = Ridge(alpha=10.0)
            r.fit(X_train[ok], y[ok])
            pred[:, j] = r.predict(X_pred)
            n_ok += 1
        except Exception:
            continue
        if (j+1) % 200 == 0:
            log(f"  {j+1}/{len(drugs)}: {n_ok} fits OK")

    pred_df = pd.DataFrame(pred, index=mps, columns=drugs).T
    pred_df = pred_df.dropna(how="all")
    log(f"  predictions: {pred_df.shape}")

    pred_df.to_csv(FIG3/"drug_sensitivity_per_mp.csv")

    # Selectivity: std across MPs (high = more state-specific)
    pred_df["std_across_mp"] = pred_df[mps].std(axis=1)
    pred_df["min_AUC"] = pred_df[mps].min(axis=1)
    pred_df["selective_for"] = pred_df[mps].idxmin(axis=1)  # lowest AUC = most sensitive
    rows = []
    for mp in mps:
        sub = pred_df[pred_df["selective_for"] == mp].sort_values(mp, ascending=True).head(20)
        for drug, vals in sub.iterrows():
            rows.append({"MP": mp, "drug": drug, "predicted_AUC": float(vals[mp]),
                          "std_across_mp": float(vals["std_across_mp"])})
    pd.DataFrame(rows).to_csv(FIG3/"drug_selectivity_summary.csv", index=False)
    log(f"  drug_selectivity_summary.csv ({len(rows)} rows)")

    with open(FIG3/"drug_sensitivity_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 20c — Drug sensitivity per MP via DepMap PRISM ridge\n\n")
        f.write(f"- Cell lines used: {len(common)}\n")
        f.write(f"- Drugs with successful ridge: {n_ok}/{len(drugs)}\n")
        f.write(f"- Common genes: {len(common_genes)}\n")
        f.write("Lower predicted AUC = more sensitive. selective_for = MP "
                "with lowest predicted AUC.\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
