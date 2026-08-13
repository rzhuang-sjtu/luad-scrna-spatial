# Where each input comes from

Two kinds of file are read by the scripts but not produced by them. Both are
listed here so that a reader who follows the pipeline does not stop at a missing
file and wonder whether a script is absent.

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
| `${DATA_ROOT}/High-resolution/neutrophil_final.h5ad` | Salcher et al. 2022 NSCLC neutrophil reference | the atlas released with that paper |
| `${PROJECT_ROOT}/data/reference/gavish2023_MPs.csv` | the 41 pan-cancer meta-programmes | supplementary table of Gavish et al. 2023 |
| Hallmark `.gmt` | MSigDB Hallmark gene sets | MSigDB, after accepting its licence |
| Geneformer V2-104M weights | the foundation model | Hugging Face, `ctheodoris/Geneformer` at the commit pinned in `environments/*-pip-freeze.txt` |

## Derived tables deposited with the paper

The following are read by plotting or table-building scripts but have no
producer in this repository. They were written during the analysis by code that
is not part of the released set, so they are distributed as derived data in the
Zenodo deposit that accompanies this repository (DOI in `CITATION.cff`).
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
