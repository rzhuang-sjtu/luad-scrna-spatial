"""Step 2.3: ST cohort × reference gene overlap, plus c2l-style gene filter."""
import os
from pathlib import Path
import scanpy as sc, numpy as np, pandas as pd

REF = "${DATA_ROOT}/ST/results/step02_reference/unified_reference.h5ad"
ST  = "${DATA_ROOT}/ST/results/step01_qc/luad_tumor_sections.h5ad"
OUT = Path("${DATA_ROOT}/ST/results/step02_reference")

ref = sc.read_h5ad(REF, backed="r")
st  = sc.read_h5ad(ST,  backed="r")
print(f"reference: {ref.shape}; ST cohort: {st.shape}")

ref_genes = set(ref.var_names)
st_genes  = set(st.var_names)
overlap   = ref_genes & st_genes
print(f"reference unique genes: {len(ref_genes)}")
print(f"ST cohort unique genes: {len(st_genes)}")
print(f"intersection: {len(overlap)}")

# c2l-style filter: cells expressing each gene in reference
# Re-load reference fully for filtering check
ref.file.close()
ref_full = sc.read_h5ad(REF)
import scipy.sparse as sp
X = ref_full.X
if sp.issparse(X):
    n_cells_per_gene = (X > 0).sum(axis=0).A1
    nonz_mean = np.asarray(X.mean(axis=0)).ravel()  # mean across cells (incl zeros)
    nonzero_only_mean = np.divide(np.asarray(X.sum(axis=0)).ravel(), n_cells_per_gene, where=n_cells_per_gene>0)
else:
    n_cells_per_gene = (X > 0).sum(axis=0)
    nonzero_only_mean = X.sum(axis=0) / np.where(n_cells_per_gene>0, n_cells_per_gene, 1)

n_cells_total = ref_full.n_obs
pct = n_cells_per_gene / n_cells_total

# c2l defaults: cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
cell_count_cutoff = 5
cell_percentage_cutoff2 = 0.03
nonz_mean_cutoff = 1.12

mask_count = n_cells_per_gene >= cell_count_cutoff
mask_pct   = pct >= cell_percentage_cutoff2
mask_nonz  = nonzero_only_mean >= nonz_mean_cutoff
# c2l rule: keep if (mask_count AND (mask_pct OR mask_nonz))
c2l_keep = mask_count & (mask_pct | mask_nonz)
ref_kept_genes = set(ref_full.var_names[c2l_keep])
print(f"\nc2l filter (count>={cell_count_cutoff} & (pct>=3% | nonz_mean>=1.12)): keep {c2l_keep.sum()} / {ref_full.n_vars}")

final_overlap = ref_kept_genes & st_genes
print(f"after c2l filter, ref ∩ ST = {len(final_overlap)}")
print(f"\n=> use {len(final_overlap)} genes for cell2location")

# Save lists
out_df = pd.DataFrame({
    "gene": list(ref_full.var_names[c2l_keep]),
    "in_ST": [g in st_genes for g in ref_full.var_names[c2l_keep]],
})
out_df.to_csv(OUT / "c2l_genes_kept.csv", index=False)

# Spot-check key Fig7 genes survival
KEY = ["OSM","OSMR","LIFR","IL6ST","IL1A","IL1B","IL1R1","SPP1","CD44","MMP9",
       "CXCL8","CXCL1","CXCL2","PLAU","PLAUR","TGFB1","HIF1A","JAK1","STAT3",
       "NFKBIA","ATF3","FOSB","JUN","JUNB","C1QC","MARCO","FOLR2","S100A8","S100A9"]
for g in KEY:
    in_ref = g in ref_genes
    in_st  = g in st_genes
    in_c2l = g in ref_kept_genes
    in_final = g in final_overlap
    print(f"   {g:8s}  ref={in_ref}  ST={in_st}  c2l_kept={in_c2l}  final={in_final}")
