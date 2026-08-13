"""Smoke test: run only the neu transition (smallest) end-to-end."""
import os, time, sys
from pathlib import Path

os.environ.setdefault("TRANSFORMERS_NO_ADVISORY_WARNINGS", "1")
os.environ.setdefault("HF_DATASETS_VERBOSITY", "error")

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats

ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer")
MODEL_DIR = "${PROJECT_ROOT}/data/external/geneformer/Geneformer-V2-104M"
T = "neu_osm_priming_to_low"
CS = {"state_key": "cell_state", "start_state": "sender", "goal_state": "receiver", "alt_states": []}

tok_dir = ROOT / "tokenized" / T / f"{T}.dataset"
pert_dir = ROOT / "perturb" / T
state_dir = pert_dir / "state_embs"
pert_dir.mkdir(parents=True, exist_ok=True)
state_dir.mkdir(parents=True, exist_ok=True)

t0 = time.time()
print("[A] state embs")
emb = EmbExtractor(emb_mode="cls", cell_emb_style="mean_pool",
                   max_ncells=2000, summary_stat="exact_mean",
                   forward_batch_size=16, nproc=4, model_version="V2")
state_embs = emb.get_state_embs(CS, MODEL_DIR, str(tok_dir), str(state_dir), T)
print("  embs:", {k: tuple(v.shape) for k, v in state_embs.items()})

print("[B] perturb (max 200 sender cells for smoke)")
isp = InSilicoPerturber(perturb_type="delete", genes_to_perturb="all",
                        emb_mode="cls", cell_emb_style="mean_pool",
                        cell_states_to_model=CS, state_embs_dict=state_embs,
                        max_ncells=200, forward_batch_size=16, nproc=4,
                        model_version="V2", clear_mem_ncells=100)
isp.perturb_data(MODEL_DIR, str(tok_dir), str(pert_dir), T)

print("[C] stats")
stats = InSilicoPerturberStats(mode="goal_state_shift", cell_states_to_model=CS)
stats.get_stats(str(pert_dir), None, str(pert_dir), f"{T}_stats")

print(f"\nELAPSED {(time.time()-t0)/60:.1f} min")
print("outputs:")
for p in pert_dir.rglob("*"):
    if p.is_file():
        print(" ", p, p.stat().st_size)
