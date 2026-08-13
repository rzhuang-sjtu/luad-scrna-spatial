#!/usr/bin/env python
"""
D1b — rewrite cell-type labels for decontX input to fix the wrong annotation file used in D1.

Background: D1 used celltype_coarse from `luad_merged_annotated.h5ad` as decontX z.
      That annotation is incorrect — of the 240,884 cells labelled Fibroblast,
      58.0% express PTPRC, 45.8% express LYZ, whereas only 10.5% express COL1A1;
      in `luad_copykat.h5ad`, Fibroblast COL1A1 positivity is 74.2% and PTPRC only 5.0%.
      The two annotations agree on only 67.96% of 853,469 cells, and the composition is clearly implausible
      (bad annotation: Fibroblast 28%, B cells 0.15%).

decontX z directly determines the background model; wrong labels yield wrong estimates and must be re-run.
Matrices and barcodes are unchanged; only labels.tsv is rewritten — avoids re-exporting 131 sample matrices.

Usage: python D1b_relabel_decontx.py
Then re-run D2_run_decontx.R.
"""
import os, glob, shutil
import numpy as np
import pandas as pd
import h5py
import anndata as ad

IN = "${PROJECT_ROOT}/results/decontx/input"
CK = "${PROJECT_ROOT}/data/processed/luad_copykat.h5ad"
NEU = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad"

print("Loading correct labels (luad_copykat.h5ad)...", flush=True)
with h5py.File(CK, "r") as h:
    ix = h["obs"].attrs.get("_index", "index")
    cid = np.array([x.decode() for x in h["obs"][ix][:]])
    g = h["obs"]["celltype_coarse"]
    cats = np.array([x.decode() for x in g["categories"][:]])
    cd = g["codes"][:]
    ct = np.where(cd >= 0, cats[np.clip(cd, 0, None)], "Unknown")
lab = pd.Series(ct, index=cid)
print(f"{len(lab):,} cells; composition {lab.value_counts().to_dict()}", flush=True)

print("Loading neutrophil barcodes ...", flush=True)
neu_bc = set(ad.read_h5ad(NEU, backed="r").obs_names)
print(f"{len(neu_bc):,}", flush=True)

samples = sorted(os.path.basename(d) for d in glob.glob(f"{IN}/*") if os.path.isdir(d))
print(f"samples: {len(samples)}", flush=True)

rows = []
for i, s in enumerate(samples, 1):
    d = f"{IN}/{s}"
    bc = pd.read_csv(f"{d}/barcodes.tsv", header=None)[0].astype(str).values
    old = pd.read_csv(f"{d}/labels.tsv", header=None)[0].astype(str).values
    assert len(bc) == len(old), f"{s}: barcode and label row counts differ"

    new = lab.reindex(bc).fillna("Unknown").values.astype(object)
    isneu = np.array([b in neu_bc for b in bc])
    new[isneu] = "Neutrophil"

    # Back up original labels once for later comparison (do not overwrite if already present)
    if not os.path.exists(f"{d}/labels_merged_annotated.tsv"):
        shutil.copy(f"{d}/labels.tsv", f"{d}/labels_merged_annotated.tsv")
    pd.Series(new).to_csv(f"{d}/labels.tsv", index=False, header=False)

    rows.append(dict(sample=s, n_cells=len(bc),
                     n_neutrophil=int(isneu.sum()),
                     n_unknown=int((new == "Unknown").sum()),
                     n_types_old=len(set(old)), n_types_new=len(set(new)),
                     agree=float((old == new).mean())))
    if i % 25 == 0:
        print(f"  {i}/{len(samples)}", flush=True)

R = pd.DataFrame(rows)
R.to_csv("${PROJECT_ROOT}/results/decontx/D1b_relabel_summary.csv", index=False)
print(f"\nOld–new label agreement: median {R.agree.median():.3f}"
      f"（{R.agree.min():.3f} ~ {R.agree.max():.3f}）", flush=True)
print(f"Cells without a matched label: total {R.n_unknown.sum():,}", flush=True)
print(f"n types old median {R.n_types_old.median():.0f} → new median {R.n_types_new.median():.0f}",
      flush=True)
print("Original labels backed up as labels_merged_annotated.tsv; re-run D2_run_decontx.R now", flush=True)
