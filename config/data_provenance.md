# Where each input comes from

Two kinds of file are read by the scripts but not produced by them. The tables
below list the ones a reader will hit, so that following the pipeline does not
stop at a missing file and leave you wondering whether a script is absent. If you
find an input that is read but not listed here, it is an omission — please open
an issue.

## External data you download yourself

The repository contains no data. Each accession below is public; the paper cites
them individually.

| Expected location | What it is | Where to get it |
|---|---|---|
| `${DATA_ROOT}/<GSE accession>/` | the seven scRNA-seq cohorts | GEO, accessions in the README |
| `${DATA_ROOT}/ST/E-MTAB-13530/` | discovery Visium sections | ArrayExpress E-MTAB-13530 |
| `${DATA_ROOT}/ST/Okamura 2024/Visium_FF_LUAD_No_*.tar.gz` | validation Visium sections | processed data released with Takano et al. 2024 (DOI 10.1038/s41467-024-54671-7); the raw records JGAS000613/JGAS000677 are controlled-access and are not used |
| `${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv`, `TCGA_LUAD_clinical.csv` | TCGA-LUAD expression and clinical | UCSC Xena, TCGA-LUAD open-access layer |
| `${WORK_ROOT}/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz`, `..._all_thresholded.by_genes.gz` | TCGA-LUAD GISTIC2 copy number | UCSC Xena |
| `${WORK_ROOT}/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz` | TCGA MC3 non-silent mutation matrix | UCSC Xena, MC3 public release |
| `${DATA_ROOT}/High-resolution/neutrophil_final.h5ad` | Salcher et al. 2022 NSCLC neutrophil reference | CELLxGENE <https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287> or Zenodo <https://doi.org/10.5281/zenodo.6411867>; rename/copy to `neutrophil_final.h5ad` |
| `${DATA_ROOT}/GSE223503/41591_2025_3530_MOESM3_ESM.xlsx` | sample metadata for GSE223503, read by `data_prep/qc/qc_GSE223503.py` | supplementary file MOESM3 of the Nature Medicine article that released GSE223503; download from the publisher and place under `${DATA_ROOT}/GSE223503/` |
| `${PROJECT_ROOT}/data/reference/collectri_symbols.tsv` | CollecTRI regulons behind the Fig. 3 TF-activity panel, read by `analysis/trajectory_tf/20a_tf_activity.py` | export from decoupleR/OmniPath: `decoupler.get_collectri(organism="human", split_complexes=False)` written to TSV |
| `${PROJECT_ROOT}/data/reference/wilkerson_LAD_centroids.csv` | Wilkerson LAD expression-subtype centroids, read by `analysis/cnmf_meta_programs/08b_wilkerson_subtype.py` | supplementary material of Wilkerson et al., *Clin Cancer Res* 2012 (the 506-gene LAD predictor) |
| `${PROJECT_ROOT}/data/reference/gavish2023_MPs.csv` | the 41 pan-cancer meta-programmes | supplementary tables of Gavish et al. 2023, Nature <https://doi.org/10.1038/s41586-023-06130-4>, also 3CA <https://www.weizmann.ac.il/sites/3CA>; export/rename to `gavish2023_MPs.csv` |
| `${DATA_ROOT}/depmap/24Q2/` — `Model.csv`, `CRISPRGeneEffect.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | DepMap Public 24Q2 CRISPR gene effect, cell-line model table and expression | depmap.org data page, 24Q2 release; read by six scripts under `analysis/perturbation_cascade/` and `revision/` |
| `${DATA_ROOT}/GSE207422/` — `GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz`, `..._metadata.xlsx` | neoadjuvant chemo-immunotherapy cohort | GEO GSE207422 supplementary files |
| `${DATA_ROOT}/GSE126044/GSE126044_counts.txt.gz` | anti-PD-1 responder cohort, raw counts | GEO GSE126044 supplementary file |
| `${DATA_ROOT}/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz` | anti-PD-1 durable-benefit cohort | GEO GSE135222 supplementary file |
| `${DATA_ROOT}/GSE68465/GSE68465_series_matrix.txt.gz` | Director's Challenge microarray validation cohort | GEO GSE68465 series matrix |
| `${DATA_ROOT}/GSE31210/GSE31210_family.soft.gz` | Okayama microarray validation cohort | GEO GSE31210 SOFT family file |
| `${PROJECT_ROOT}/data/gmt/MSigDB_Hallmark_2020.gmt` | Hallmark gene sets | Enrichr library `MSigDB_Hallmark_2020`, or MSigDB `h.all.*.symbols.gmt` renamed to that filename; no script downloads it |
| Geneformer V2-104M weights | the foundation model | Hugging Face, `ctheodoris/Geneformer` at the commit pinned in `environments/*-pip-freeze.txt` |

## Derived tables deposited with the paper

The following are read by plotting or table-building scripts but have no
producer in this repository. They were written during the analysis by code that
is not part of the released set, so they are distributed as derived data in the
Zenodo deposit that accompanies this repository: https://doi.org/10.5281/zenodo.21912264
Download the deposit and place them where the consuming script expects them.

| File | Read by |
|---|---|
| `dpt_root_sensitivity/composition_by_ventile.csv` | `plotting/revision_panels/P2_fig23s1_panels.R` (Fig. 2K) |
| `tf_activity_mp_zscore.csv` | `plotting/fig3/render_fig3_v3.R` and three supplementary-table scripts |
| `salcher_neu_disease.csv` | `revision/ambient_and_neutrophil/F4_luad_only_reference.py` |
| `pseudotime_mp_score_curves.csv`, `pseudotime_umap_winsorized.csv.gz` | `plotting/fig2/figure2_main.R` |
| `hallmark_gsea_top20_per_mp.csv` | `plotting/figS5/figure_s5_part1.R`, `part2.R` |
| `fig5e_ligand_mean_raw.csv`, `fig5e_rename_map.csv` | `plotting/fig5/plot_fig5_v2.R` |
| `8O_km_*_stats.csv`, `8P_km_*_stats.csv` | `plotting/fig8/panels_R_part1.R` and the supplementary-table builders |

`s1_lisi_scores.csv` and `s1_lisi_per_cell.csv.gz` used to be in this list; the
script that writes them is now included, at
`revision/integration_qc/rerun_lisi.py`.

## Intermediates the released code does produce

| Path | Written by | Read by |
|---|---|---|
| `${PROJECT_ROOT}/data/processed/*.h5ad` | `data_prep/atlas_build/` | most of `analysis/` |
| `${PROJECT_ROOT}/results/` | `analysis/`, `revision/` | `plotting/` |
| `${WORK_ROOT}/luad_figures/` | `analysis/`, `revision/` | `plotting/` |
| `${DATA_ROOT}/<GSE>/*_cleaned.h5ad` | `data_prep/qc/qc_GSE*.py` | `data_prep/atlas_build/01–03` |

One naming mismatch is worth knowing about: the QC scripts write
`<GSE>_cleaned.h5ad` while `data_prep/atlas_build/01_inspect_metadata.py` reads
`<GSE>_clean.h5ad` from `${WORK_ROOT}/数据清洗/`. The cleaned files were copied
and renamed into that directory by hand between the two stages. Either rename on
copy, or adjust the paths at the top of the atlas_build scripts.

## A third staging directory

`data_prep/qc/dataset_inclusion_analysis.py` — the script README points to as the
justification for excluding GSE149655, GSE223503 and GSE308103 — reads the QC
outputs from a third location, `${DATA_ROOT}/cleaned/`, under the names listed at
the top of that file (`<GSE>_clean.h5ad` or `<GSE>_LUAD_clean.h5ad`). Nothing in
this repository writes that directory: copy the QC outputs there by hand. The
script also rewrites two of those files in place.

## The figure-data tree

Analysis scripts write their plot inputs under `${PROJECT_ROOT}/results/`
(`fig5_plot_data/`, `fig8_plot_data/v2_500/`, and so on). The plotting scripts and
the supplementary-table builders read them from a second tree,
`${WORK_ROOT}/luad_figures/<figure>/`, under different directory names — for
example `results/fig5_plot_data/` versus `luad_figures/fig5/data/`. Nothing in
this repository copies between the two: that step was done by hand.

To run the plotting or table-building scripts, copy each producer's output into
the matching `luad_figures/` directory first. `analysis/supplementary_tables/`
uses `safe_read`, which prints `[SKIP] missing <path>` and omits the sheet rather
than failing, so a workbook built without the copy step will be silently short of
sheets — check the console output against the sheet list in the manuscript.
