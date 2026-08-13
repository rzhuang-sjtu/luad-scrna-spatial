"""Prepare Geneformer V2 input h5ads for the 3 LUAD state transitions (Fig 8A).

Output structure:
  ~/luad/results/fig8_geneformer/inputs/{transition}/data.h5ad
  with var.ensembl_id, obs.n_counts, obs.cell_state in {'sender','receiver'}.

Transitions:
  macro_spp1_to_c1qc : myeloid h5ad, dominant_macro = argmax(sc_Macro_*)
  mal_mp3_to_mp1     : malignant h5ad, dominant_MP
  neu_osm_priming_to_low : neu h5ad, neu_subtype
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

GFM_DIR = Path("${PROJECT_ROOT}/data/external/geneformer/geneformer")
OUT_DIR = Path("${PROJECT_ROOT}/results/fig8_geneformer/inputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

with open(GFM_DIR / "gene_name_id_dict_gc104M.pkl", "rb") as f:
    NAME2ID = pickle.load(f)
with open(GFM_DIR / "token_dictionary_gc104M.pkl", "rb") as f:
    TOK = pickle.load(f)


def map_ensembl(adata):
    """Add var.ensembl_id and drop genes without mapping or not in vocab."""
    var = adata.var.copy()
    var["ensembl_id"] = [NAME2ID.get(s, "") for s in adata.var_names]
    keep = [eid in TOK and eid != "" for eid in var["ensembl_id"]]
    print(f"  [map_ensembl] kept {sum(keep)}/{len(keep)} genes (in V2 vocab)")
    a = adata[:, keep].copy()
    a.var = var.loc[keep].copy()
    a.var_names = a.var["ensembl_id"].values
    a.var.index.name = None
    return a


def get_counts_layer(adata):
    """Return counts as sparse CSR; prefer 'counts' layer, fallback to X if integer."""
    if "counts" in adata.layers:
        X = adata.layers["counts"]
    elif "count" in adata.layers:
        X = adata.layers["count"]
    else:
        X = adata.X
        if hasattr(X, "toarray"):
            sample = X[:50].toarray()
        else:
            sample = X[:50]
        if not np.allclose(sample, sample.astype(int)):
            raise ValueError("X is not integer counts and no counts layer found")
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    return X.astype(np.float32).tocsr()


def write_transition(name, adata_full, sender_mask, receiver_mask):
    """Slice + write input h5ad for a transition."""
    print(f"\n=== {name} ===")
    print(f"  sender={sender_mask.sum()} receiver={receiver_mask.sum()}")
    keep = sender_mask | receiver_mask
    a = adata_full[keep].copy()
    cs = np.where(sender_mask[keep], "sender", "receiver")
    a.obs["cell_state"] = pd.Categorical(cs, categories=["sender", "receiver"])
    # counts + n_counts
    X = get_counts_layer(a)
    a.X = X
    for k in list(a.layers.keys()):
        del a.layers[k]
    if a.raw is not None:
        a.raw = None
    a.obs["n_counts"] = np.asarray(X.sum(axis=1)).flatten().astype(np.float32)
    # ensembl mapping
    a = map_ensembl(a)
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out / "data.h5ad")
    print(f"  wrote {out/'data.h5ad'}: shape={a.shape}, n_counts median={np.median(a.obs.n_counts):.0f}")


# === transition 1: Macro_SPP1 → Macro_C1QC ===
print("Loading myeloid h5ad...")
mye = sc.read_h5ad("${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad")
sc_cols = [
    "sc_Macro_general", "sc_Macro_C1QC", "sc_Macro_SPP1", "sc_Macro_FCN1",
    "sc_Macro_FOLR2", "sc_Macro_MARCO", "sc_Macro_prolif",
]
M = mye.obs[sc_cols].astype(float).values
labels = np.array([c.replace("sc_", "") for c in sc_cols])
dom = labels[np.argmax(M, axis=1)]
macro_mask = mye.obs["myeloid_subtype"].astype(str).str.startswith("Macro_").values
dom = np.where(macro_mask, dom, "")
mye.obs["dominant_macro"] = pd.Categorical(dom)
sender = (mye.obs["dominant_macro"] == "Macro_SPP1").values
receiver = (mye.obs["dominant_macro"] == "Macro_C1QC").values
write_transition("macro_spp1_to_c1qc", mye, sender, receiver)

# === transition 2: Mal_MP3 → Mal_MP1 ===
print("\nLoading malignant h5ad...")
mal = sc.read_h5ad("${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad")
sender = (mal.obs["dominant_MP"] == "MP3").values
receiver = (mal.obs["dominant_MP"] == "MP1").values
write_transition("mal_mp3_to_mp1", mal, sender, receiver)

# === transition 3: Neu_OSM_priming → Neu_OSM_low ===
print("\nLoading neutrophil h5ad...")
neu = sc.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad")
sender = (neu.obs["neu_subtype"] == "Neu_OSM_priming").values
receiver = (neu.obs["neu_subtype"] == "Neu_OSM_low").values
write_transition("neu_osm_priming_to_low", neu, sender, receiver)

print("\nDONE.")
