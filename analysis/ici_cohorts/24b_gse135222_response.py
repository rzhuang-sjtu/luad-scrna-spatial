"""Step 24b: GSE135222 anti-PD-1/L1 NSCLC (n=27) — MP score vs DCB.

Response definition (DCB convention):
  R (DCB)   = PFS ≥ 180 days OR censored (no progression observed)
  NR (NDB)  = PFS < 180 days AND event = 1 (progressed)
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

GSE = "gse135222"
EXPR_FILE = Path("${DATA_ROOT}/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")
DCB_THRESHOLD_DAYS = 180


def ensg_to_symbol(ensg_ids: list[str]) -> dict[str, str]:
    """Convert ENSG (versioned) → official gene symbol via mygene."""
    import mygene
    mg = mygene.MyGeneInfo()
    stripped = [e.split(".")[0] for e in ensg_ids]
    res = mg.querymany(stripped, scopes="ensembl.gene", fields="symbol",
                        species="human", returnall=False, verbose=False)
    mapping = {}
    for r in res:
        if r.get("notfound"):
            continue
        mapping[r["query"]] = r.get("symbol", "")
    out = {}
    for orig, stripped_id in zip(ensg_ids, stripped):
        sym = mapping.get(stripped_id, "")
        if sym:
            out[orig] = sym
    return out


def main():
    C.log(f"=== {GSE.upper()} anti-PD-1/L1 NSCLC ===")
    expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)
    C.log(f"raw shape (ENSG × samples): {expr.shape}")
    expr.index.name = "ensg"

    # Map ENSG → Symbol
    cache = Path.home() / "luad/data/external/gse135222_ensg2symbol.csv"
    if cache.exists():
        m = pd.read_csv(cache).set_index("ensg")["symbol"].to_dict()
        C.log(f"loaded cached ENSG→Symbol: {len(m)}")
    else:
        m = ensg_to_symbol(expr.index.tolist())
        pd.DataFrame({"ensg": list(m.keys()),
                      "symbol": list(m.values())}).to_csv(cache, index=False)
        C.log(f"fetched ENSG→Symbol: {len(m)}, cached")
    expr = expr.loc[[i for i in expr.index if i in m]].copy()
    expr.index = [m[i] for i in expr.index]
    expr.index.name = "gene"
    expr = expr.groupby(level=0).max().astype("float32")
    C.log(f"gene-level shape: {expr.shape}")

    # Per the GEO summary, file is RNA-seq expression (likely RPKM/TPM).
    # Header of file shows decimals like 16.86, 34.47 — looks like TPM/RPKM.
    # Apply log2(x+1).
    expr = np.log2(expr + 1).astype("float32")

    # Clinical
    gse_obj = C.fetch_geo_soft("GSE135222")
    meta = C.gsm_characteristics(gse_obj)
    meta["sample_id"] = meta["title"].str.replace(" ", "")  # "NSCLC 990" → "NSCLC990"
    meta = meta.set_index("sample_id")
    meta["pfs_event"] = pd.to_numeric(meta["progression-free_survival_(pfs)"], errors="coerce").astype(int)
    meta["pfs_days"] = pd.to_numeric(meta["pfs.time"], errors="coerce").astype(float)
    common = [s for s in expr.columns if s in meta.index]
    C.log(f"matched samples: {len(common)}")
    expr = expr[common]; meta = meta.loc[common]

    # Build DCB groups
    is_R = ((meta["pfs_event"] == 0) | (meta["pfs_days"] >= DCB_THRESHOLD_DAYS))
    is_NR = (meta["pfs_event"] == 1) & (meta["pfs_days"] < DCB_THRESHOLD_DAYS)
    meta["response_group"] = "ambig"
    meta.loc[is_R.values, "response_group"] = "R"
    meta.loc[is_NR.values, "response_group"] = "NR"
    C.log(f"groups: {meta['response_group'].value_counts().to_dict()}")

    # ssGSEA
    gene_sets, overlaps = C.build_gene_sets(expr.index)
    for k, n in overlaps.items():
        C.log(f"  {k}: {n} present")
    C.log("running ssGSEA")
    scores = C.run_ssgsea(expr, gene_sets)
    C.log(f"  scores shape: {scores.shape}")

    score_df = scores.join(meta[["response_group", "pfs_event", "pfs_days",
                                  "gender", "age"]], how="inner")
    score_df.index.name = "Sample"
    # Drop ambiguous for the comparison
    score_df_for_cmp = score_df[score_df["response_group"].isin(["R","NR"])].copy()
    C.log(f"R={int((score_df_for_cmp['response_group']=='R').sum())}, "
          f"NR={int((score_df_for_cmp['response_group']=='NR').sum())}")

    comp = C.write_outputs(GSE, score_df_for_cmp, group_col="response_group",
                            pos_label="R", neg_label="NR")
    # Re-write full mp_scores including ambiguous (for record)
    score_df.to_csv(C.OUT / f"{GSE}_mp_scores.csv")
    C.log("\n=== Comparison (R=DCB ≥180d positive) ===")
    print(comp.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
