"""Step 10d: GO BP + KEGG enrichment per myeloid subtype DEG list.

Uses step10_myeloid_markers.csv (per-subtype top DEGs). For each subtype,
pick genes with logFC>0.5 and pval_adj<0.05 (fallback: top 100 by score).
Run gseapy.enrichr with GO_Biological_Process_2023 and KEGG_2021_Human.
Top 10 terms per gene set per subtype.

Outputs:
  ~/luad/results/step10_myeloid_go_enrichment.csv
  ${WORK_ROOT}/luad_figures/fig4/myeloid_go_enrichment.csv
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import pandas as pd
import gseapy as gp

MARKERS = Path.home() / "luad/results/step10_myeloid_markers.csv"
RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig4")
GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human"]


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {MARKERS}")
    markers = pd.read_csv(MARKERS)
    log(f"  DEG table: {markers.shape}, subtypes: {markers['subtype'].unique()}")

    out_rows = []
    for sub, d in markers.groupby("subtype"):
        # Preferred: logFC>0.5 & pval_adj<0.05
        tight = d[(d["logFC"] > 0.5) & (d["pval_adj"] < 0.05)].copy()
        if len(tight) >= 10:
            gene_list = tight["gene"].dropna().astype(str).tolist()
            source = "logFC>0.5 & padj<0.05"
        else:
            gene_list = d.sort_values("score", ascending=False).head(100)["gene"].astype(str).tolist()
            source = "top 100 by score"
        log(f"  {sub}: {len(gene_list)} genes ({source})")
        if len(gene_list) < 5:
            continue
        try:
            enr = gp.enrichr(gene_list=gene_list, gene_sets=GENE_SETS,
                             organism="human", outdir=None, no_plot=True)
            for _, row in enr.results.iterrows():
                out_rows.append({
                    "subtype": sub,
                    "gene_set": row["Gene_set"],
                    "term": row["Term"],
                    "overlap": row["Overlap"],
                    "p_value": row["P-value"],
                    "adj_p_value": row["Adjusted P-value"],
                    "odds_ratio": row.get("Odds Ratio", None),
                    "combined_score": row.get("Combined Score", None),
                    "genes": row["Genes"],
                    "source_filter": source,
                })
        except Exception as e:
            log(f"    enrichr failed for {sub}: {e}")

    df = pd.DataFrame(out_rows)
    log(f"  all enrichr rows: {len(df)}")
    # keep all but also produce top-10 per subtype×gene_set by combined_score
    if len(df) > 0:
        top10 = (df.sort_values("combined_score", ascending=False, na_position="last")
                    .groupby(["subtype", "gene_set"]).head(10).reset_index(drop=True))
        df.to_csv(RES/"step10_myeloid_go_enrichment_full.csv", index=False)
        top10.to_csv(RES/"step10_myeloid_go_enrichment.csv", index=False)
        top10.to_csv(FIG/"myeloid_go_enrichment.csv", index=False)
        log("  top-3 terms per subtype×set:")
        log(top10.groupby(["subtype","gene_set"]).head(3)[
            ["subtype","gene_set","term","adj_p_value","combined_score"]
        ].round(4).to_string(index=False))
    else:
        log("  no enrichment results to write")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
