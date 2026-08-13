"""R1 — rebuild the Neu_OSM_priming -> Neu_OSM_low Geneformer input on the FULL gene matrix.

Bug being fixed
---------------
01_prepare_geneformer_inputs.py loaded `luad_neutrophil_own_annotated.h5ad`, which had
already been subset to 3,000 highly-variable genes by the scVI/scANVI step. The macrophage
and malignant transitions were built from full-gene objects (9,422 / 9,603 genes), so only
the neutrophil arm was truncated:

    transition                  genes in input      mean tokens/cell
    macro_spp1_to_c1qc                   9,422               1,385.5
    mal_mp3_to_mp1                       9,603               2,028.3
    neu_osm_priming_to_low               2,982                  94.6   <-- truncated

Measured on the same 1,880 cells, the full matrix (9,698 genes; 9,480 in the V2 vocabulary)
gives a mean of 380.9 tokens/cell (median 354, min 173, max 1,106) — a 4.0x recovery.

Counts come from luad_neutrophil_own_raw.h5ad (X is raw integer counts, same cell order as
the annotated object); labels come from the annotated object.

Writes to inputs_v2/ so the published inputs/ tree is left untouched.
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

GFM_DIR = Path("${PROJECT_ROOT}/data/external/geneformer/geneformer")
OUT_DIR = Path("${PROJECT_ROOT}/results/fig8_geneformer/inputs_v2")
OUT_DIR.mkdir(parents=True, exist_ok=True)

RAW = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad"
ANN = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad"
SENDER, RECEIVER = "Neu_OSM_priming", "Neu_OSM_low"

# split-half seed: two disjoint random halves of the sender pool, written as separate
# inputs so the split is random rather than length-stratified (InSilicoPerturber sorts
# by token length internally, so cell_inds_to_perturb cannot give a random split).
SPLIT_SEED = 20260804

NAME2ID = pickle.load(open(GFM_DIR / "gene_name_id_dict_gc104M.pkl", "rb"))
TOK = pickle.load(open(GFM_DIR / "token_dictionary_gc104M.pkl", "rb"))


def map_ensembl(adata):
    var = adata.var.copy()
    var["ensembl_id"] = [NAME2ID.get(s, "") for s in adata.var_names]
    keep = [eid in TOK and eid != "" for eid in var["ensembl_id"]]
    print(f"  [map_ensembl] kept {sum(keep)}/{len(keep)} genes (in V2 vocab)")
    a = adata[:, keep].copy()
    a.var = var.loc[keep].copy()
    a.var_names = a.var["ensembl_id"].values
    a.var.index.name = None
    return a


def write_input(name, adata_full, sender_mask, receiver_mask):
    print(f"\n=== {name} ===")
    print(f"  sender={int(sender_mask.sum())} receiver={int(receiver_mask.sum())}")
    keep = sender_mask | receiver_mask
    a = adata_full[keep].copy()
    a.obs["cell_state"] = pd.Categorical(
        np.where(sender_mask[keep], "sender", "receiver"), categories=["sender", "receiver"]
    )
    X = a.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.astype(np.float32).tocsr()
    assert np.allclose(X[:200].toarray(), np.round(X[:200].toarray())), "X is not integer counts"
    a.X = X
    for k in list(a.layers.keys()):
        del a.layers[k]
    a.raw = None
    a.obs["n_counts"] = np.asarray(X.sum(axis=1)).flatten().astype(np.float32)
    a = map_ensembl(a)
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out / "data.h5ad")
    nz = np.asarray((a.X > 0).sum(axis=1)).ravel()
    print(f"  wrote {out/'data.h5ad'}  shape={a.shape}")
    print(f"  detected genes/cell: mean={nz.mean():.1f} median={np.median(nz):.0f} "
          f"min={nz.min()} max={nz.max()}")


print("loading neutrophil objects ...")
raw = sc.read_h5ad(RAW)
ann = sc.read_h5ad(ANN)
assert (raw.obs_names == ann.obs_names).all(), "cell order mismatch between raw and annotated"
raw.obs["neu_subtype"] = ann.obs["neu_subtype"].values
print(f"  raw {raw.shape}  (published input used {ann.shape[1]} HVGs)")

sender = (raw.obs["neu_subtype"] == SENDER).values
receiver = (raw.obs["neu_subtype"] == RECEIVER).values

# full input (used for the 500-cell rerun and the saturation curve)
write_input("neu_osm_priming_to_low", raw, sender, receiver)

# two disjoint random halves of the sender pool for split-half reproducibility
rng = np.random.default_rng(SPLIT_SEED)
idx = np.where(sender)[0]
perm = rng.permutation(idx)
h1, h2 = perm[: len(perm) // 2], perm[len(perm) // 2:]
for tag, half in [("neu_osm_priming_to_low_half1", h1), ("neu_osm_priming_to_low_half2", h2)]:
    m = np.zeros(len(sender), bool)
    m[half] = True
    write_input(tag, raw, m, receiver)

print("\nDONE. Next: R2_tokenize.py")
