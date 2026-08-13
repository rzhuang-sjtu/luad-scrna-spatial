"""G1 — Export per-patient inputs for inferCNV (independent CNV cross-check).

Concern addressed: non-epithelial aneuploid cells were treated as false positives and forced non-malignant, with no independent CNV cross-validation
Methodologically problematic. This pipeline rechecks CopyKAT malignant calls with independent inferCNV.

Reference-cell strategy (closest to CopyKAT paired mode; also the official inferCNV recommendation):
  Observation cells = all epithelial cells from that patient
  Reference cells = the patient's own T/NK and B cells (genomically normal), capped at 500
  Skip patients with <50 reference cells; report coverage honestly in the results

Export format (standard inferCNV trio):
  counts.mtx      raw counts, genes x cells
  genes.tsv       gene names (aligned to gene_order)
  cells.tsv       cell barcodes
  annot.tsv       two columns: barcode \t group ("epithelial" or "reference")
  Also export gene_order.tsv at the top level (gene chromosome start end) for inferCNV ordering
"""
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import anndata as ad, os, subprocess

OUT = "${PROJECT_ROOT}/results/infercnv/input"
os.makedirs(OUT, exist_ok=True)
MAX_REF, MIN_REF, MIN_EPI = 500, 50, 50

print("Export gene locus table (using CopyKAT built-in hg38 annotation)...", flush=True)
GO = f"{OUT}/../gene_order.tsv"
os.makedirs(os.path.dirname(GO), exist_ok=True)
subprocess.run(["Rscript", "-e", f'''
suppressMessages(library(copykat))
a <- get("full.anno", envir=asNamespace("copykat"))
a <- a[!is.na(a$hgnc_symbol) & a$hgnc_symbol != "", ]
a <- a[!duplicated(a$hgnc_symbol), ]
a$chr <- paste0("chr", a$chromosome_name)
write.table(a[, c("hgnc_symbol","chr","start_position","end_position")],
            "{GO}", sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("gene_order rows:", nrow(a), "\\n")
'''], check=True)

print("Loading merged object ...", flush=True)
a = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad")
ct = a.obs["celltype_coarse"].astype(str).values
pt = (a.obs["dataset"].astype(str) + "__" + a.obs["patient_id"].astype(str)).values
X = sp.csr_matrix(a.X)
genes = np.array(a.var_names)
go = pd.read_csv(GO, sep="\t", header=None, names=["g", "chr", "s", "e"])
keep_g = np.isin(genes, go.g.values)
print(f"{a.shape}; genes with locus annotation {keep_g.sum()}", flush=True)
X = X[:, keep_g]; genes = genes[keep_g]

rng = np.random.default_rng(42)
rows, skipped = [], []
pts = sorted(set(pt))
for i, p in enumerate(pts):
    m = pt == p
    epi = np.where(m & (ct == "Epithelial"))[0]
    ref = np.where(m & np.isin(ct, ["T_NK", "B"]))[0]
    if len(epi) < MIN_EPI or len(ref) < MIN_REF:
        skipped.append((p, len(epi), len(ref))); continue
    if len(ref) > MAX_REF:
        ref = rng.choice(ref, MAX_REF, replace=False)
    idx = np.concatenate([epi, ref])
    lab = np.array(["epithelial"] * len(epi) + ["reference"] * len(ref))
    d = f"{OUT}/{p}"; os.makedirs(d, exist_ok=True)
    sio.mmwrite(f"{d}/counts.mtx", X[idx].T.tocoo())
    pd.Series(genes).to_csv(f"{d}/genes.tsv", index=False, header=False)
    bc = a.obs_names[idx]
    pd.Series(bc).to_csv(f"{d}/cells.tsv", index=False, header=False)
    pd.DataFrame({"bc": bc, "grp": lab}).to_csv(f"{d}/annot.tsv", sep="\t",
                                                index=False, header=False)
    rows.append(dict(patient=p, n_epi=len(epi), n_ref=len(ref), n_total=len(idx)))
    if (i + 1) % 20 == 0:
        print(f"{i+1}/{len(pts)} processed, exported {len(rows)}", flush=True)

R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/../patients.csv", index=False)
print(f"\nExported {len(R)} patients; skipped {len(skipped)} (epithelial<{MIN_EPI} or reference<{MIN_REF})", flush=True)
print(f"Total epithelial cells {R.n_epi.sum()}, total reference cells {R.n_ref.sum()}", flush=True)
print(f"Epithelial per patient median {R.n_epi.median():.0f} (max {R.n_epi.max()}),"
      f"reference median {R.n_ref.median():.0f}", flush=True)
if skipped:
    print(f"First 5 skipped: {skipped[:5]}", flush=True)
sz = subprocess.run(["du", "-sh", OUT], capture_output=True, text=True).stdout.split()[0]
print(f"Export directory size {sz}", flush=True)
