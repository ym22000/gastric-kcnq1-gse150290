# Gastric_KCNQ1

This folder contains the gastric single-cell analysis kept for `dataset 1` only, based on `GSE150290`.

The question was to test whether `KCNQ1`, `KCNE2`, and `KCNE3` are associated with a broad intestinal-versus-diffuse tumor axis in gastric cancer.

## Folder structure

- `bibliography/`
  Source paper and supplementary methods.
- `pipeline_markers_only/`
  Main analysis branch.
- `pipeline_markers_plus_cnv/`
  Same dataset with an additional CNV-based refinement step.
- `cache/`
  Processed GEO matrices used in the main workflow.
- `raw_geo/`
  Raw GEO files kept locally for traceability and CNV reruns.
- `deps/`
  Local dependency files used during setup.

## Pipeline summary

The workflow was reproduced as closely as possible from the paper using the parameters reported in the article and supplementary methods.

The analysis was done in two branches:

### 1. `pipeline_markers_only`

- build the global non-immune atlas from the processed GEO matrices
- identify the broad tumor compartment
- re-cluster tumor cells with paper-like Seurat settings
- organize tumor cells along an intestinal-versus-diffuse axis
- project `KCNQ1`, `KCNE2`, and `KCNE3`

### 2. `pipeline_markers_plus_cnv`

- keep the same global logic
- add an inferCNV-based validation layer
- refine the tumor compartment with marker and CNV support
- re-evaluate `KCNQ1`, `KCNE2`, and `KCNE3` in the refined tumor map

## Paper fidelity

This repository follows the published workflow as closely as possible from the material available locally.

It is not exactly the same pipeline as the original study because the authors' full code and detailed internal annotations were not available. The reconstruction therefore relies on the paper, the supplementary methods, and the released matrices.

## Main result

Across the two branches, the interpretation stayed consistent:

- `KCNQ1` is more intestinal-oriented
- `KCNE3` follows the same direction
- `KCNE2` is weaker

## Main files

Figures:

- `pipeline_markers_only/outputs/tumor_poles/figures/tumor_all_global_intestinal_diffuse_annotated.pdf`
- `pipeline_markers_only/outputs/kcnq1/figures/KCNQ1_featureplot_global_and_tumor.pdf`
- `pipeline_markers_only/outputs/kcne2_kcne3/figures/KCNE2_featureplot_global_and_tumor.pdf`
- `pipeline_markers_only/outputs/kcne2_kcne3/figures/KCNE3_featureplot_global_and_tumor.pdf`
- `pipeline_markers_plus_cnv/outputs/global_cnv_screen/figures/global_tsne_tumor_markers_and_cnv_like_genes.pdf`
- `pipeline_markers_plus_cnv/outputs/refined_kcnq1_kcne2_kcne3/figures/refined_tumor_KCNQ1_KCNE2_KCNE3_multiplot.pdf`

Scripts:

- `pipeline_markers_only/scripts/run_global_nonimmune_tumor_extraction.R`
- `pipeline_markers_only/scripts/run_all_global_tumor_cells_exact_paper_pipeline.R`
- `pipeline_markers_only/scripts/annotate_broad_tumor_tsne_intestinal_diffuse.R`
- `pipeline_markers_plus_cnv/scripts/run_infercnv_on_existing_global_tsne.R`
- `pipeline_markers_plus_cnv/scripts/refine_tumor_compartment_markers_plus_cnv.R`
- `pipeline_markers_plus_cnv/scripts/plot_refined_kcnq1_kcne2_kcne3_and_spearman.R`
