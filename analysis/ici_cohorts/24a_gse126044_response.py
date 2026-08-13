"""Step 24a: GSE126044 anti-PD-1 NSCLC (n=16) — MP score vs R/NR."""
from __future__ import annotations
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(os.path.expanduser("~/luad/scripts")))
from pathlib import Path
import numpy as np
import pandas as pd

from importlib import import_module
C = import_module("24_common")

GSE = "gse126044"
COUNTS = Path("${DATA_ROOT}/GSE126044/GSE126044_counts.txt.gz")


def main():
    C.log(f"=== {GSE.upper()} anti-PD-1 NSCLC ===")
    # 1. counts
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
    counts.index.name = "gene"
    if not counts.index.is_unique:
        counts = counts.groupby(level=0).max()
    C.log(f"counts shape: {counts.shape}; samples={list(counts.columns)}")

    # 2. clinical from SOFT
    gse_obj = C.fetch_geo_soft("GSE126044")
    meta = C.gsm_characteristics(gse_obj)
    meta["sample_id"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    meta = meta.set_index("sample_id")
    common_samples = [s for s in counts.columns if s in meta.index]
    C.log(f"matched samples: {len(common_samples)}")
    counts = counts[common_samples]
    meta = meta.loc[common_samples]

    # 3. log2(CPM+1) (gseapy ssgsea is rank-based but consistent w/ TCGA pipeline)
    libsize = counts.sum(axis=0)
    cpm = (counts.div(libsize, axis=1) * 1e6).astype("float32")
    expr = np.log2(cpm + 1).astype("float32")
    C.log(f"log2(CPM+1) shape: {expr.shape}")

    # 4. gene sets
    gene_sets, overlaps = C.build_gene_sets(expr.index)
    for k, n in overlaps.items():
        C.log(f"  {k}: {n} present")

    # 5. ssGSEA
    C.log("running ssGSEA")
    scores = C.run_ssgsea(expr, gene_sets)
    C.log(f"  scores shape: {scores.shape}")

    # 6. join clinical
    score_df = scores.join(meta[["patient_response", "sample"]], how="inner")
    score_df["response_group"] = score_df["patient_response"].map(
        {"responder": "R", "non-responder": "NR"}).fillna(score_df["patient_response"])
    score_df.index.name = "Sample"

    # 7. write
    comp = C.write_outputs(GSE, score_df, group_col="response_group",
                            pos_label="R", neg_label="NR")
    C.log("\n=== Comparison (R=responder positive) ===")
    print(comp.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
