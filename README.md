# Analysis code — LUAD single-cell, spatial and in-silico perturbation study

Code for *Deep-learning integration of single-cell and spatial transcriptomics
reveals a neutrophil-associated mesenchymal niche and a candidate translocon
dependency in lung adenocarcinoma* (Cancer Immunology, Immunotherapy).

This repository contains **code only**. Every dataset the released workflow
reads is publicly accessible; accessions are listed below and cited individually
in the paper. Raw records for the validation spatial cohort (JGAS000613,
JGAS000677) are controlled-access: the workflow does not read them and they are
not redistributed here — it reads the openly available processed Visium data.
The Zenodo deposit <https://doi.org/10.5281/zenodo.21912264> holds two things:
an archive of this repository, and a `derived_tables/` bundle with the tables
that the plotting scripts read but this repository does not regenerate (file
list in `config/data_provenance.md`). The same concept DOI is used for both
because they are versioned together.

## Layout

```
data_prep/    Build the h5ad objects and input tables from raw downloads.
                atlas_build/         merge the seven cohorts, integrate
                qc/                  per-dataset quality control
                copykat/             malignant-cell calling
                cnmf/                inputs for consensus NMF
                neutrophil/          reference-based label transfer inputs
                cellchat_inputs/     objects for communication analysis
                spatial/             Visium QC and deconvolution reference
                perturbation/        Geneformer inputs and tokenisation
                plot_data_exports/   per-figure source tables
analysis/     Computation.
                annotation/          cell-type assignment
                cnmf_meta_programs/  malignant meta-programmes
                trajectory_tf/       pseudotime and regulon activity
                myeloid_neutrophil/  myeloid and neutrophil subsets
                neutrophil_emt_link/ neutrophil states versus malignant
                                     programmes, EMT ligands
                tnk/                 T and NK subsets
                cellchat/            communication networks
                liana_nets/          ligand-receptor aggregation
                spatial/             MP scoring, PROGENy, COMMOT, MISTy
                tcga_survival/       Cox models and Kaplan-Meier
                external_validation/ independent expression cohorts
                ici_cohorts/         immunotherapy response
                perturbation_cascade/ in-silico knock-out and the
                                     multi-endpoint validation cascade
                supplementary_tables/ build the supplementary workbooks
plotting/     Rendering, one directory per published figure: fig1-fig8 and
              figS1-S7, figS9, figS11. Figure numbers are kept here on purpose:
              they are how a reader locates the code behind a panel. Figure S10
              is rendered by plotting/fig7/stepD_fig_s10_plot.R; figures S8 and
              S12 have no plotting code of their own.
                revision_panels/     the panels added during revision
                                     (Fig. 2K, 3K, 8N-8P, S1E/S1F, S2D, S13)
                shared/              theme and palette used across figures
revision/     Analyses added during revision, grouped by topic:
                geneformer_stability/    rank reproducibility, null cascade
                copy_number_and_cnv/     7p11.2 confounding, inferCNV
                survival/                Cox diagnostics, canonical re-fits
                ambient_and_neutrophil/  decontX, neutrophil subtypes
                spatial_and_trajectory/  ROI thresholds, depth adjustment,
                                         pseudotime root sensitivity
                integration_qc/          kBET / LISI
                mp3_sec61g_robustness/   model-specification checks
              The panels these analyses produce are rendered by
              plotting/revision_panels/.
config/       Path configuration (config/paths.md) and where every input
              comes from (config/data_provenance.md).
environments/ Conda environments, pip freezes, R package list.
```

## Where to start

To reproduce a specific panel, open the matching directory under `plotting/`;
each script names the source table it reads, and those tables are written by
`data_prep/plot_data_exports/`. To follow an analysis end to end, start from
`data_prep/atlas_build/` and work through `analysis/` in the order the layout
above lists.

## Running the code

1. Create the environments in `environments/` (see `INSTALL.md` there).
2. Set the four path roots described in `config/paths.md`.
3. Download the datasets listed below to `${DATA_ROOT}`.
4. Within each folder, scripts run in numerical order.

Scripts were written for a single workstation and are not packaged as a
pipeline; each is standalone and reads and writes files under the configured
roots.

## Datasets

### Single-cell RNA-seq (LUAD atlas)

| Accession | Use |
|---|---|
| GSE123902 | atlas merge |
| GSE131907 | atlas merge |
| GSE143423 | atlas merge (brain metastases) |
| GSE148071 | atlas merge |
| GSE164789 | atlas merge (precursor lesions) |
| GSE189357 | atlas merge |
| GSE253013 | atlas merge |
| GSE149655, GSE223503, GSE308103 | passed QC but not merged into the atlas; the inclusion decision is in `data_prep/qc/dataset_inclusion_analysis.py` |
| Salcher et al. 2022 NSCLC atlas — CELLxGENE <https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287> or Zenodo <https://doi.org/10.5281/zenodo.6411867>; place the neutrophil reference at `${DATA_ROOT}/High-resolution/neutrophil_final.h5ad`, renaming if needed | scANVI reference for neutrophil label transfer |

### Spatial transcriptomics

| Accession | Use |
|---|---|
| E-MTAB-13530 | discovery cohort, 12 sections |
| JGAS000613 / JGAS000677 (Takano et al.) | validation cohort, 8 sections (controlled access; processed Visium openly available) |

### Bulk expression and survival

| Accession | Use |
|---|---|
| TCGA-LUAD (GDC) | survival, copy number (GISTIC2), mutations (MC3) |
| GSE68465 | external validation |
| GSE31210 | external validation |

### Immunotherapy response

| Accession | Use |
|---|---|
| GSE207422 | neoadjuvant anti-PD-1 |
| GSE126044 | anti-PD-1 |
| GSE135222 | anti-PD-1 |
| GSE93157 | screened, not used — too few LUAD cases (`analysis/ici_cohorts/24c_gse93157_skip.py`) |
| GSE14814, GSE42127 | adjuvant chemotherapy cohorts, screened for the composite score |
| GSE243013 | screened; metadata only |

### Other

| Resource | Use |
|---|---|
| DepMap Public 24Q2 | CRISPR gene effect, copy number, cell-line expression |
| Geneformer V2-104M (`ctheodoris/Geneformer`, commit ad8f66d) | in-silico perturbation |

## Citation

See `CITATION.cff`. Paper: *Cancer Immunology, Immunotherapy*. Code archive: https://doi.org/10.5281/zenodo.21912264

## Contact

Code: Ruizhe Huang, <rhuang6@outlook.com>  
Correspondence: Siyu Chen, <siyu.chen@shsmu.edu.cn>

## Network activity

Some scripts reach the internet. All of it is either downloading public data or
querying a public annotation service; nothing uploads data, and no expression
matrix, sample table or file ever leaves the machine.

| What | Where | Sent | Received |
|---|---|---|---|
| Public data download | NCBI GEO, Figshare (DepMap), MSigDB | an accession | the dataset |
| CellTypist model | the CellTypist model store | nothing | model weights |
| Enrichment testing | Enrichr | a list of gene symbols | enrichment results |
| Gene ID mapping | MyGene | a list of Ensembl IDs | matching symbols |

MSigDB asks you to accept its licence before downloading. No script fetches the
Hallmark sets for you: download them yourself (Enrichr library
`MSigDB_Hallmark_2020`, or MSigDB `h.all.*.symbols.gmt` renamed to match) and
place the file at `${PROJECT_ROOT}/data/gmt/MSigDB_Hallmark_2020.gmt` — every
consumer reads that exact path and filename.

## License

MIT — see `LICENSE`. The MIT terms cover the code in this repository only.
Runtime dependencies are installed by the user and keep their own licences;
several are GPL (for example NMF, CellChat, CopyKAT, edgeR, limma, survival).
No dependency source is redistributed here, so calling them from these scripts
does not place the repository under GPL.
