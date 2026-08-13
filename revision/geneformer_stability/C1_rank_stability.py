"""C1 — Geneformer rank stability (split-half / bootstrap / saturation); no GPU required.

Motivation: a perturbation ranking is only usable if it is reproducible.
This measures that directly — split-half and bootstrap agreement at 500
sender cells, and a saturation curve out to >=1,000.

Rationale
----
InSilicoPerturber pickle outputs are stored **per cell** (the
cell_embs_{i} index is the cell index). Perturbation results for one cell do not depend on which other cells are in the batch
(model in eval mode; no batchnorm/dropout; padding masked); the target-state centroid is also computed independently.
Passing any cell-subset pickles to InSilicoPerturberStats is therefore equivalent to re-running on that subset alone.

Verified: recomputing on all 500 cells yields Shift_to_goal_end and FDR identical to the published table
bit-for-bit (max difference 0.0); all 2,278 genes match.

Thus split-half, bootstrap, and saturation curves with n<=500 need no new GPU time.
Only points with n>500 require a real run (see run_neu_local.sh).

The two curves measure different things; report both:
  consistency(n) — top-200 Jaccard between two **non-overlapping** n-cell subsets.
                     Direct measure of whether n cells suffice; can be extrapolated.
  agreement(n)   — top-200 Jaccard between an n-cell subset and the full 500-cell result.
                     Rises as n approaches 500; only reflects subset fraction and cannot alone argue convergence.
"""
import os, re, glob, json, shutil, tempfile, time, argparse
import numpy as np
import pandas as pd
from scipy import stats as sps
from geneformer import InSilicoPerturberStats

ROOT = "${PROJECT_ROOT}/results/fig8_geneformer"
OUT = "${PROJECT_ROOT}/results/fig8_geneformer/rank_stability"
CELL_STATES = {"state_key": "cell_state", "start_state": "sender",
               "goal_state": "receiver", "alt_states": []}
FDR, N_MIN, TOPN = 0.05, 5, 200


def cell_files(src):
    """{cell_index: [pickle paths]}"""
    d = {}
    for f in glob.glob(os.path.join(src, "*_raw.pickle")):
        m = re.search(r"cell_embs_(\d+)batch", os.path.basename(f))
        if m:
            d.setdefault(int(m.group(1)), []).append(f)
    return d


def stats_for(cells, cf, dup_tag=None):
    """Run InSilicoPerturberStats on the given cell subset; return the filtered ranking."""
    tmp = tempfile.mkdtemp(prefix="gfstab_")
    try:
        for k, c in enumerate(cells):
            for f in cf[c]:
                b = os.path.basename(f)
                if dup_tag is not None:      # bootstrap: same cell may be drawn more than once; rename to keep unique
                    b = re.sub(r"cell_embs_(\d+)batch", f"cell_embs_{9000000+k}batch", b)
                dst = os.path.join(tmp, b)
                if not os.path.exists(dst):
                    os.symlink(f, dst)
        st = InSilicoPerturberStats(mode="goal_state_shift", genes_perturbed="all", combos=0,
                                    anchor_gene=None, cell_states_to_model=CELL_STATES,
                                    pickle_suffix="_raw.pickle")
        st.get_stats(input_data_directory=tmp, null_dist_data_directory=None,
                     output_directory=tmp, output_prefix="s")
        d = pd.read_csv(os.path.join(tmp, "s.csv"), index_col=0)
        d = d[(d.Sig == 1) & (d.Shift_to_goal_end > 0) &
              (d.Goal_end_FDR < FDR) & (d.N_Detections >= N_MIN)]
        return d.sort_values("Shift_to_goal_end", ascending=False).reset_index(drop=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a or b) else np.nan


def rank_corr(da, db):
    m = da[["Gene_name", "Shift_to_goal_end"]].merge(
        db[["Gene_name", "Shift_to_goal_end"]], on="Gene_name", suffixes=("_a", "_b"))
    if len(m) < 20:
        return np.nan, len(m)
    return sps.spearmanr(m.Shift_to_goal_end_a, m.Shift_to_goal_end_b)[0], len(m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturb-dir", default=f"{ROOT}/perturb_500")
    ap.add_argument("--transitions", nargs="*", default=[
        "macro_spp1_to_c1qc", "mal_mp3_to_mp1", "neu_osm_priming_to_low"])
    ap.add_argument("--n-splithalf", type=int, default=10)
    ap.add_argument("--n-boot", type=int, default=10)
    ap.add_argument("--sat-reps", type=int, default=3)
    ap.add_argument("--tag", default="published500")
    # Sampling points for the consistency curve (two non-overlapping subsets need 2n <= n_cells); increase when more cells are available
    ap.add_argument("--sat-ns", type=int, nargs="*",
                    default=[25, 50, 100, 150, 200, 250])
    ap.add_argument("--agree-ns", type=int, nargs="*",
                    default=[50, 100, 200, 300, 400, 450])
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(20260804)
    rows_sh, rows_sat, rows_boot, rows_ref = [], [], [], []

    for t in a.transitions:
        src = os.path.join(a.perturb_dir, t)
        cf = cell_files(src)
        cells = np.array(sorted(cf))
        print(f"\n=== {t}: {len(cells)} cells ===", flush=True)
        t0 = time.time()
        ref = stats_for(cells, cf)
        ref_top = list(ref.Gene_name.head(TOPN))
        rows_ref.append({"transition": t, "n_cells": len(cells), "n_pass_filter": len(ref),
                         "seconds": round(time.time() - t0, 1)})
        print(f"  reference: {len(ref)} genes pass filter ({time.time()-t0:.0f}s)", flush=True)

        # --- split-half ---
        for r in range(a.n_splithalf):
            p = rng.permutation(cells); h = len(cells) // 2
            A, B = stats_for(p[:h], cf), stats_for(p[h:2 * h], cf)
            rho, nov = rank_corr(A, B)
            rows_sh.append({"transition": t, "rep": r, "n_per_half": h,
                            "jaccard_top200": jac(A.Gene_name.head(TOPN), B.Gene_name.head(TOPN)),
                            "spearman_shift": rho, "n_common_genes": nov,
                            "n_pass_A": len(A), "n_pass_B": len(B)})
            print(f"  split-half {r}: J={rows_sh[-1]['jaccard_top200']:.3f} rho={rho:.3f}", flush=True)

        # --- saturation ---
        for n in a.sat_ns:
            if 2 * n > len(cells): continue
            for r in range(a.sat_reps):
                p = rng.permutation(cells)
                A, B = stats_for(p[:n], cf), stats_for(p[n:2 * n], cf)
                rows_sat.append({"transition": t, "n": n, "rep": r, "kind": "disjoint_consistency",
                                 "jaccard_top200": jac(A.Gene_name.head(TOPN), B.Gene_name.head(TOPN)),
                                 "spearman_shift": rank_corr(A, B)[0]})
            print(f"  consistency n={n}: "
                  f"J={np.mean([x['jaccard_top200'] for x in rows_sat if x['n']==n and x['transition']==t]):.3f}",
                  flush=True)
        for n in a.agree_ns:
            if n >= len(cells): continue
            for r in range(a.sat_reps):
                s = rng.choice(cells, n, replace=False)
                A = stats_for(s, cf)
                rows_sat.append({"transition": t, "n": n, "rep": r, "kind": "agreement_with_full",
                                 "jaccard_top200": jac(A.Gene_name.head(TOPN), ref_top),
                                 "spearman_shift": rank_corr(A, ref)[0]})

        # --- bootstrap ---
        for r in range(a.n_boot):
            s = rng.choice(cells, len(cells), replace=True)
            A = stats_for(s, cf, dup_tag=r)
            rows_boot.append({"transition": t, "rep": r,
                              "jaccard_top200": jac(A.Gene_name.head(TOPN), ref_top),
                              "spearman_shift": rank_corr(A, ref)[0], "n_pass": len(A)})
            print(f"  bootstrap {r}: J={rows_boot[-1]['jaccard_top200']:.3f}", flush=True)

    for name, rows in [("reference", rows_ref), ("splithalf", rows_sh),
                       ("saturation", rows_sat), ("bootstrap", rows_boot)]:
        pd.DataFrame(rows).to_csv(f"{OUT}/{a.tag}_{name}.csv", index=False)

    print("\n" + "=" * 70)
    print(" RANK STABILITY SUMMARY (top-200 Jaccard)")
    print("=" * 70)
    sh = pd.DataFrame(rows_sh); sat = pd.DataFrame(rows_sat); bo = pd.DataFrame(rows_boot)
    for t in a.transitions:
        s = sh[sh.transition == t]
        b = bo[bo.transition == t]
        print(f"\n{t}")
        print(f"  split-half (250 vs 250) : J = {s.jaccard_top200.mean():.3f} "
              f"± {s.jaccard_top200.std():.3f}   Spearman = {s.spearman_shift.mean():.3f}")
        print(f"  bootstrap  (vs full 500): J = {b.jaccard_top200.mean():.3f} ± {b.jaccard_top200.std():.3f}")
        c = sat[(sat.transition == t) & (sat.kind == "disjoint_consistency")]
        if len(c):
            print("  consistency between disjoint subsets:")
            for n, g in c.groupby("n"):
                print(f"    n={n:>3}: J = {g.jaccard_top200.mean():.3f}")
    print(f"\nwrote CSVs to {OUT}/")


if __name__ == "__main__":
    main()
