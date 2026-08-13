"""D1 — Prepare per-sample inputs for decontX.

decontX must be run on **all cell types** from the same sample: it models each cell's expression as a mixture of
the true cell-type profile and a sample-shared background profile, and infers the background from expression that
should not be present in a given cell type. Neutrophils alone cannot identify the background.

For each sample with neutrophils, export raw counts for all cells plus cell-type labels,
with neutrophils labelled as a separate class (they were merged into Myeloid in the original celltype_coarse).

Output: results/decontx/input/{sample}/{matrix.mtx, genes.tsv, barcodes.tsv, labels.tsv}
"""
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import anndata as ad, os

OUT = "${PROJECT_ROOT}/results/decontx/input"
os.makedirs(OUT, exist_ok=True)
MIN_CELLS = 100          # decontX estimates are unstable when the sample has too few cells
MIN_NEU = 10             # require at least this many neutrophils to run

print("Loading merged object ...", flush=True)
a = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad")
print(f"  {a.shape}", flush=True)

# Neutrophil barcodes (from myeloid subset with myeloid_subtype == Neutrophil)
neu = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad", backed="r")
neu_bc = set(neu.obs_names)
print(f"neutrophils {len(neu_bc)}", flush=True)

lab = a.obs["celltype_coarse"].astype(str).values.copy()
is_neu = np.array([b in neu_bc for b in a.obs_names])
lab[is_neu] = "Neutrophil"
a.obs["dx_label"] = lab
print(f"labelled Neutrophil: {is_neu.sum()}", flush=True)

samples = pd.Series(a.obs["sample_id"].astype(str).values)
neu_per_sample = pd.Series(samples[is_neu]).value_counts()
todo = [s for s, n in neu_per_sample.items() if n >= MIN_NEU]
print(f"samples with >={MIN_NEU} neutrophils: {len(todo)}", flush=True)

X = sp.csr_matrix(a.X)
genes = np.array(a.var_names)
kept, skipped = [], []
for i, s in enumerate(todo):
    m = (samples == s).values
    if m.sum() < MIN_CELLS:
        skipped.append((s, int(m.sum()))); continue
    d = f"{OUT}/{s}"
    os.makedirs(d, exist_ok=True)
    sub = X[m]
    sio.mmwrite(f"{d}/matrix.mtx", sub.T.tocoo())      # decontX expects genes x cells
    pd.Series(genes).to_csv(f"{d}/genes.tsv", index=False, header=False)
    pd.Series(a.obs_names[m]).to_csv(f"{d}/barcodes.tsv", index=False, header=False)
    pd.Series(a.obs["dx_label"].values[m]).to_csv(f"{d}/labels.tsv", index=False, header=False)
    kept.append({"sample": s, "n_cells": int(m.sum()),
                 "n_neu": int((is_neu & m).sum()),
                 "n_types": int(pd.Series(a.obs["dx_label"].values[m]).nunique())})
    if (i + 1) % 20 == 0:
        print(f"exported {i+1}/{len(todo)}", flush=True)

K = pd.DataFrame(kept)
K.to_csv("${PROJECT_ROOT}/results/decontx/samples.csv", index=False)
print(f"\nExported {len(K)} samples, skipped {len(skipped)} (cells <{MIN_CELLS})", flush=True)
print(f"total cells {K.n_cells.sum()}, of which neutrophils {K.n_neu.sum()}", flush=True)
print(f"cells per sample median {K.n_cells.median():.0f}, neutrophils median {K.n_neu.median():.0f}", flush=True)
print(f"cell types per sample median {K.n_types.median():.0f}", flush=True)
