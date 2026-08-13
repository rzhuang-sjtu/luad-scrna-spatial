"""Run Geneformer V2 in-silico KO for the 3 LUAD transitions and rank genes.

For each transition:
  A. EmbExtractor.get_state_embs → state_embs_dict {sender: emb, receiver: emb}
  B. InSilicoPerturber with cell_states_to_model + state_embs_dict
     → raw pickle of per-gene embedding-shift toward receiver
  C. InSilicoPerturberStats → ranked CSV of genes that shift sender → receiver
"""
import os
import sys
import time
import pickle
from pathlib import Path

# silence transformers/datasets noise
os.environ.setdefault("TRANSFORMERS_NO_ADVISORY_WARNINGS", "1")
os.environ.setdefault("HF_DATASETS_VERBOSITY", "error")

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats

ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer")
MODEL_DIR = "${PROJECT_ROOT}/data/external/geneformer/Geneformer-V2-104M"

TRANSITIONS = ["macro_spp1_to_c1qc", "mal_mp3_to_mp1", "neu_osm_priming_to_low"]
CELL_STATES = {
    "state_key": "cell_state",
    "start_state": "sender",
    "goal_state": "receiver",
    "alt_states": [],
}

MAX_SENDER_CELLS = 500   # subsample sender cells per transition (controls runtime)
FWD_BS = 16               # forward batch size; V2 input_size=4096; RTX 3080 20GB (32 caused cudaErrorNotReady on mal long-seq)

if __name__ == "__main__":
 for t in TRANSITIONS:
    t0 = time.time()
    print(f"\n========== transition: {t} ==========")
    tok_dir = ROOT / "tokenized" / t / f"{t}.dataset"
    pert_dir = ROOT / "perturb_500" / t
    pert_dir.mkdir(parents=True, exist_ok=True)
    state_dir = pert_dir / "state_embs"
    state_dir.mkdir(parents=True, exist_ok=True)

    # ----- A. State embeddings -----
    print(f"  [A] extracting state embs ({t})")
    emb = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cls",            # V2 uses CLS token
        cell_emb_style="mean_pool",
        max_ncells=2000,           # cap cells per state for speed; uses random subsample
        emb_layer=-1,
        summary_stat="exact_mean",
        forward_batch_size=FWD_BS,
        nproc=4,
        model_version="V2",
    )
    state_embs = emb.get_state_embs(
        cell_states_to_model=CELL_STATES,
        model_directory=MODEL_DIR,
        input_data_file=str(tok_dir),
        output_directory=str(state_dir),
        output_prefix=t,
        output_torch_embs=True,
    )
    print(f"    states: {list(state_embs.keys())}; emb dim: {state_embs['sender'].shape}")

    # ----- B. In-silico perturb (delete each gene in sender cells, measure shift) -----
    print(f"  [B] running InSilicoPerturber ({t})")
    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb="all",
        combos=0,
        anchor_gene=None,
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cls",
        cell_emb_style="mean_pool",
        cell_states_to_model=CELL_STATES,
        state_embs_dict=state_embs,
        max_ncells=MAX_SENDER_CELLS,
        cell_inds_to_perturb="all",
        emb_layer=-1,
        forward_batch_size=FWD_BS,
        nproc=4,
        model_version="V2",
        clear_mem_ncells=200,
    )
    isp.perturb_data(
        model_directory=MODEL_DIR,
        input_data_file=str(tok_dir),
        output_directory=str(pert_dir),
        output_prefix=t,
    )

    # ----- C. Stats / ranking -----
    print(f"  [C] computing perturbation stats ({t})")
    stats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed="all",
        combos=0,
        anchor_gene=None,
        cell_states_to_model=CELL_STATES,
        pickle_suffix="_raw.pickle",
    )
    stats.get_stats(
        input_data_directory=str(pert_dir),
        null_dist_data_directory=None,
        output_directory=str(pert_dir),
        output_prefix=f"{t}_stats",
    )

    elapsed = (time.time() - t0) / 60
    print(f"  [done] {t} elapsed: {elapsed:.1f} min")

print("\nALL TRANSITIONS DONE.")
