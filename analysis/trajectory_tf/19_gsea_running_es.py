"""Step 19: export GSEA running enrichment score (RES) curves for Fig S5.

Reconstructs per-MP gene rankings from cNMF pool (mean z-score across member
GEPs), reruns gseapy.prerank, and extracts the running enrichment score
(RES) vector + hit positions per Hallmark term.

Outputs to ${WORK_ROOT}/luad_figures/fig2/:
  - gsea_running_es.csv.gz       long format: MP, Term, rank, running_es, is_hit
  - gsea_hit_positions.csv       MP, Term, hit_rank (stem markers only)
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd

POOL = Path.home()/"luad/data/cnmf_output/gep_pool_zscore.csv"
ASSIGN = Path.home()/"luad/results/step6_gep_mp_assignment.csv"
GMT = Path.home()/"luad/data/gmt/MSigDB_Hallmark_2020.gmt"
OUT = Path("${WORK_ROOT}/luad_figures/fig2")
PERM_N = 10  # minimal permutations; FDR is taken from existing step7 output


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log("loading cNMF pool + MP assignments")
    pool = pd.read_csv(POOL, index_col=0)
    log(f"  pool: {pool.shape}")
    assign = pd.read_csv(ASSIGN)
    mp_col = "MP" if "MP" in assign.columns else \
             [c for c in assign.columns if c != "gep_id"][0]
    log(f"  assign: {assign.shape}, MP col={mp_col}")

    # Per-MP ranking: mean z-score across member GEPs, descending
    log("building per-MP rankings")
    rankings = {}
    for mp, sub in assign.groupby(mp_col):
        if mp == "MP5": continue
        members = [g for g in sub["gep_id"].tolist() if g in pool.columns]
        if len(members) < 5: continue
        mean_z = pool[members].mean(axis=1).sort_values(ascending=False)
        rankings[str(mp)] = mean_z
        log(f"  {mp}: {len(members)} GEPs, ranked {len(mean_z)} genes")

    log(f"running gp.prerank (perm_num={PERM_N}) for {len(rankings)} MPs")
    import gseapy as gp
    es_rows = []
    hit_rows = []
    for mp, ranked_series in rankings.items():
        t1 = time.time()
        rnk = pd.DataFrame({"gene": ranked_series.index,
                              "rank": ranked_series.values})
        pre = gp.prerank(
            rnk=rnk, gene_sets=str(GMT), outdir=None,
            min_size=5, max_size=1000,
            permutation_num=PERM_N, seed=0, threads=4, verbose=False,
        )
        # gseapy 1.x: pre.results is dict[term] → dict with RES, hits, etc.
        n_genes = len(ranked_series)
        for term, info in pre.results.items():
            res_array = np.asarray(info.get("RES", []))
            if res_array.size == 0:
                continue
            hits_idx = list(info.get("hits", info.get("hits_indices", [])))
            # Append running ES (entire curve)
            for r, val in enumerate(res_array, start=1):
                es_rows.append((mp, term, r, float(val)))
            # Hit positions for stem markers
            for h in hits_idx:
                hit_rows.append((mp, term, int(h) + 1))
        log(f"  {mp}: {len(pre.results)} terms in {time.time()-t1:.1f}s")

    log(f"writing output files (rows: ES={len(es_rows)}, hits={len(hit_rows)})")
    es_df = pd.DataFrame(es_rows, columns=["MP", "Term", "rank", "running_es"])
    es_df.to_csv(OUT/"gsea_running_es.csv.gz", index=False, compression="gzip")
    log(f"  gsea_running_es.csv.gz {es_df.shape}  size={(OUT/'gsea_running_es.csv.gz').stat().st_size/1e6:.1f} MB")

    hit_df = pd.DataFrame(hit_rows, columns=["MP", "Term", "hit_rank"])
    hit_df.to_csv(OUT/"gsea_hit_positions.csv", index=False)
    log(f"  gsea_hit_positions.csv {hit_df.shape}")

    # Sanity: print ES curve summary for one example term
    ex_mp, ex_term = "MP1", "TNF-alpha Signaling via NF-kB"
    sub = es_df[(es_df["MP"]==ex_mp) & (es_df["Term"]==ex_term)]
    if len(sub):
        log(f"\n{ex_mp} / {ex_term}: ES curve length={len(sub)}, "
            f"max={sub['running_es'].max():.3f} at rank "
            f"{sub.loc[sub['running_es'].idxmax(), 'rank']}, "
            f"min={sub['running_es'].min():.3f}")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
