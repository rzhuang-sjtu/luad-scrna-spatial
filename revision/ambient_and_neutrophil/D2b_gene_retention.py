#!/usr/bin/env python
"""
D2b — summarise retention of key genes in neutrophils after decontX.

Supports a configurable output directory for comparing decontX runs under old vs new labels.

Three quantities per gene (neutrophils only, pooled across samples):
  det_before  detection rate before ambient RNA correction (fraction of cells with count > 0)
  det_after   detection rate after ambient RNA correction
  cnt_ret     total counts after / total counts before ambient RNA correction

Genes fall into three classes:
  target             OSM / IL1B / CXCL8 — claimed in the paper to be neutrophil-secreted
  neutrophil-intrinsic CSF3R / S100A8 — undisputed neutrophil markers; should be nearly fully retained
  ambient control    EPCAM / C1QB / APOE — should not be expressed in neutrophils; detection is mostly background and should drop sharply
The last class is an internal check that ambient RNA correction is working:
if these genes are also nearly fully retained, decontX did not estimate background and high retention of target genes is unconvincing.

Usage:
  python D2b_gene_retention.py [output_dir] [result_filename]
  Defaults: results/decontx/output and results/decontx/gene_retention.csv
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.io as sio

ROOT = "${PROJECT_ROOT}/results/decontx"
OUTD = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/output"
DEST = sys.argv[2] if len(sys.argv) > 2 else f"{ROOT}/gene_retention.csv"
IN = f"{ROOT}/input"

GENES = {"OSM": "target", "IL1B": "目标", "CXCL8": "目标",
         "CSF3R": "neutrophil-intrinsic", "S100A8": "中性粒自身",
         "EPCAM": "ambient control", "C1QB": "ambient对照", "APOE": "ambient对照"}

samples = sorted(os.path.basename(d) for d in glob.glob(f"{OUTD}/*")
                 if os.path.isdir(d))
print(f"samples {len(samples)}, output directory {OUTD}", flush=True)

# Per-gene accumulators: detecting cells, total counts, total cells
before_det = {g: 0 for g in GENES}
after_det = {g: 0 for g in GENES}
before_cnt = {g: 0.0 for g in GENES}
after_cnt = {g: 0.0 for g in GENES}
n_cells = 0
n_used = 0

for i, s in enumerate(samples, 1):
    fd = f"{OUTD}/{s}/decontaminated.mtx"
    fb = f"{OUTD}/{s}/neu_barcodes.tsv"
    if not (os.path.exists(fd) and os.path.exists(fb)):
        continue
    genes = pd.read_csv(f"{IN}/{s}/genes.tsv", header=None)[0].astype(str).values
    bc_all = pd.read_csv(f"{IN}/{s}/barcodes.tsv", header=None)[0].astype(str).values
    neu_bc = pd.read_csv(fb, header=None)[0].astype(str).values

    gi = {g: int(np.where(genes == g)[0][0]) for g in GENES if (genes == g).any()}
    if not gi:
        continue
    pos = pd.Index(bc_all).get_indexer(neu_bc)
    if (pos < 0).any():
        print(f"{s}: {int((pos<0).sum())} neutrophil barcodes missing from the original matrix; skip", flush=True)
        continue

    raw = sio.mmread(f"{IN}/{s}/matrix.mtx").tocsr()[:, pos]
    dec = sio.mmread(fd).tocsr()
    if dec.shape[1] != len(neu_bc):
        print(f"{s}: decontaminated matrix has {dec.shape[1]} columns != {len(neu_bc)} neutrophils; skip",
              flush=True)
        continue

    n_cells += len(neu_bc)
    n_used += 1
    for g, j in gi.items():
        rb = raw[j].toarray().ravel()
        db = dec[j].toarray().ravel()
        before_det[g] += int((rb > 0).sum()); after_det[g] += int((db > 0).sum())
        before_cnt[g] += float(rb.sum());     after_cnt[g] += float(db.sum())
    if i % 25 == 0:
        print(f"  {i}/{len(samples)}", flush=True)

rows = []
for g, cls in GENES.items():
    rows.append(dict(
        gene=g, cls=cls,
        det_before=before_det[g] / n_cells if n_cells else np.nan,
        det_after=after_det[g] / n_cells if n_cells else np.nan,
        cnt_ret=after_cnt[g] / before_cnt[g] if before_cnt[g] else np.nan))
R = pd.DataFrame(rows)
R.to_csv(DEST, index=False)
print(f"\nUsed {n_used} samples, {n_cells:,} neutrophils", flush=True)
print(R.to_string(index=False), flush=True)
print(f"\nWrote {DEST}", flush=True)
