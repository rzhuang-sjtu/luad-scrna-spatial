"""R3 — parameterised Geneformer V2 in-silico deletion runner.

Runs one (transition x n_cells) job end to end: state embeddings -> perturbation -> stats.
Designed to be launched repeatedly to build a saturation curve, and to be portable to a
rented GPU without editing the file (everything is a command-line flag).

Examples
--------
# smoke test: 2 cells, 30 genes, a couple of minutes on any GPU
python R3_perturb.py --transition neu_osm_priming_to_low --tokenized-root v2 \
       --max-ncells 2 --batch-size 32 --smoke --out-tag smoke

# real run, full-gene neutrophil, 500 sender cells
python R3_perturb.py --transition neu_osm_priming_to_low --tokenized-root v2 \
       --max-ncells 500 --batch-size 128 --out-tag n500

# saturation point on the published macrophage input
python R3_perturb.py --transition macro_spp1_to_c1qc --tokenized-root v1 \
       --max-ncells 1000 --batch-size 64 --out-tag n1000

Notes
-----
* --batch-size is the single biggest lever. The published run used 16 because the
  malignant transition has 4,096-token sequences and 20 GB could not hold more. Raise it
  on a larger card, but ramp (32 -> 64 -> 128) and watch nvidia-smi: the length
  distribution is skewed and the longest batch sets the peak.
* --min-detect restricts perturbation to genes detected in at least N sender cells.
  Genes below the threshold are dropped by the published post-hoc filter (n_det >= 5)
  anyway, so this costs no information and cuts runtime.
* Every point of a saturation curve must come from the SAME machine and environment,
  otherwise hardware/version differences are confounded with cell number.
"""
import argparse, os, json, pickle, time
from pathlib import Path

os.environ.setdefault("TRANSFORMERS_NO_ADVISORY_WARNINGS", "1")
os.environ.setdefault("HF_DATASETS_VERBOSITY", "error")

import numpy as np

ROOT = Path(os.environ.get("GF_ROOT", "${PROJECT_ROOT}/results/fig8_geneformer"))
MODEL_DIR = os.environ.get("GF_MODEL", "${PROJECT_ROOT}/data/external/geneformer/Geneformer-V2-104M")

TOK_ROOTS = {"v1": ROOT / "tokenized", "v2": ROOT / "tokenized_v2"}
INP_ROOTS = {"v1": ROOT / "inputs", "v2": ROOT / "inputs_v2"}

CELL_STATES = {"state_key": "cell_state", "start_state": "sender",
               "goal_state": "receiver", "alt_states": []}


def gene_subset(inp_h5ad, min_detect, top_n=None):
    """ENSEMBL ids of genes detected in >= min_detect sender cells.

    top_n: if given, return the top_n most frequently detected instead (smoke tests —
    an arbitrary slice can miss the few cells actually being perturbed).
    """
    import anndata as ad, scipy.sparse as sp
    a = ad.read_h5ad(inp_h5ad)
    m = (a.obs["cell_state"].astype(str) == "sender").values
    X = a.X[m]
    X = sp.csr_matrix(X) if not sp.issparse(X) else X
    det = np.asarray((X > 0).sum(axis=0)).ravel()
    if top_n is not None:
        order = np.argsort(-det)[:top_n]
        keep = [a.var_names[i] for i in order]
        print(f"  gene subset: top {len(keep)} most-detected genes "
              f"(detected in {det[order].min()}-{det[order].max()} of {m.sum()} sender cells)")
        return keep
    keep = [g for g, d in zip(a.var_names, det) if d >= min_detect]
    print(f"  gene subset: {len(keep)}/{a.n_vars} detected in >= {min_detect} sender cells")
    return keep


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--transition", required=True)
    p.add_argument("--tokenized-root", default="v1", choices=["v1", "v2"])
    p.add_argument("--max-ncells", type=int, default=500)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--emb-max-ncells", type=int, default=2000)
    p.add_argument("--min-detect", type=int, default=0,
                   help="restrict perturbation to genes detected in >=N sender cells (0 = all)")
    p.add_argument("--out-tag", required=True)
    p.add_argument("--nproc", type=int, default=1)  # HF datasets multi-process map can kill workers sporadically; default to single process
    p.add_argument("--clear-mem-ncells", type=int, default=200)
    p.add_argument("--smoke", action="store_true",
                   help="tiny end-to-end check: 30 genes only")
    args = p.parse_args()

    from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats

    t = args.transition
    tok = TOK_ROOTS[args.tokenized_root] / t / f"{t}.dataset"
    inp = INP_ROOTS[args.tokenized_root] / t / "data.h5ad"
    out = ROOT / f"perturb_{args.out_tag}" / t
    out.mkdir(parents=True, exist_ok=True)
    state_dir = out / "state_embs"; state_dir.mkdir(exist_ok=True)
    assert tok.exists(), f"missing tokenized dataset: {tok}"

    genes = "all"
    if args.smoke:
        genes = gene_subset(inp, 1, top_n=30)
        print(f"  SMOKE MODE: perturbing {len(genes)} genes only")
    elif args.min_detect > 0:
        genes = gene_subset(inp, args.min_detect)

    meta = dict(vars(args)); meta["tokenized"] = str(tok); meta["model"] = MODEL_DIR
    meta["n_genes_perturbed"] = (len(genes) if genes != "all" else "all")
    t0 = time.time()
    print(f"\n========== {t} | tag={args.out_tag} | cells={args.max_ncells} | bs={args.batch_size} ==========")

    print("  [A] state embeddings")
    emb = EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cls",
                       cell_emb_style="mean_pool", max_ncells=args.emb_max_ncells,
                       emb_layer=-1, summary_stat="exact_mean",
                       forward_batch_size=args.batch_size, nproc=args.nproc,
                       model_version="V2")
    state_embs = emb.get_state_embs(cell_states_to_model=CELL_STATES,
                                    model_directory=MODEL_DIR, input_data_file=str(tok),
                                    output_directory=str(state_dir), output_prefix=t,
                                    output_torch_embs=True)
    print(f"    states={list(state_embs.keys())} dim={tuple(state_embs['sender'].shape)}")
    t_emb = time.time() - t0

    print("  [B] in-silico deletion")
    isp = InSilicoPerturber(perturb_type="delete", perturb_rank_shift=None,
                            genes_to_perturb=genes, combos=0, anchor_gene=None,
                            model_type="Pretrained", num_classes=0, emb_mode="cls",
                            cell_emb_style="mean_pool", cell_states_to_model=CELL_STATES,
                            state_embs_dict=state_embs, max_ncells=args.max_ncells,
                            cell_inds_to_perturb="all", emb_layer=-1,
                            forward_batch_size=args.batch_size, nproc=args.nproc,
                            model_version="V2", clear_mem_ncells=args.clear_mem_ncells)
    isp.perturb_data(model_directory=MODEL_DIR, input_data_file=str(tok),
                     output_directory=str(out), output_prefix=t)
    t_pert = time.time() - t0 - t_emb

    print("  [C] stats")
    st = InSilicoPerturberStats(mode="goal_state_shift", genes_perturbed=genes, combos=0,
                                anchor_gene=None, cell_states_to_model=CELL_STATES,
                                pickle_suffix="_raw.pickle")
    st.get_stats(input_data_directory=str(out), null_dist_data_directory=None,
                 output_directory=str(out), output_prefix=f"{t}_stats")

    meta["seconds_state_embs"] = round(t_emb, 1)
    meta["seconds_perturb"] = round(t_pert, 1)
    meta["seconds_total"] = round(time.time() - t0, 1)
    (out / "run_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\n  [done] total {meta['seconds_total']/60:.1f} min "
          f"(embs {t_emb/60:.1f}, perturb {t_pert/60:.1f})")
    print(f"  outputs: {out}")


if __name__ == "__main__":
    main()
