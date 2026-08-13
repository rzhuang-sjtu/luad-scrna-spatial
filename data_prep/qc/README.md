# QC pipeline overview (Step 0, upstream of data_prep/atlas_build/01_*)

The released `qc_GSE*.py` are the notebook cells flattened into scripts. Each
performs the same six-step QC sequence on one public scRNA-seq cohort and writes
a cleaned `.h5ad` that feeds into `data_prep/atlas_build/02_merge_datasets.py`:

1. **Format unification** - convert raw input to AnnData (10x h5,
   matrix.mtx + barcodes + features, h5ad, csv).
2. **Per-dataset QC** - filter low-quality cells (min_genes, max_pct_mito,
   max_pct_ribo) and low-abundance genes (min_cells).
3. **Per-dataset doublet removal** - Scrublet (default expected_doublet_rate).
4. **Cohort merge + gene-name unification** - intersection of HGNC symbols.
5. **Harmony batch correction + integration QC** - kBET / LISI checks.
6. **Coarse clustering + residual low-quality cluster removal**.

Cleaned AnnData proceeds to integrated downstream analysis.
