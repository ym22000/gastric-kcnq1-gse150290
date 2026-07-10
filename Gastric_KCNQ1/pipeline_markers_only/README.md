# Pipeline Markers Only

This is the main analysis branch for `GSE150290`.

It is the clearest branch to present to a collaborator because it follows the biological logic of the paper while staying simple: first define a broad tumor compartment in the global atlas, then organize tumor cells along a broad intestinal-versus-diffuse axis, then project `KCNQ1`, `KCNE2`, and `KCNE3`.

## Input

Processed GEO matrices from:

- `../cache/`

Main reference documents:

- `../bibliography/article_gastrique.pdf`
- `../bibliography/supplementary_data_gastrique.pdf`

## Script order

### 1. Build the global non-immune atlas and select broad tumor clusters

Script:

- `scripts/run_global_nonimmune_tumor_extraction.R`

Main outputs:

- `outputs/global_nonimmune_atlas/figures/global_nonimmune_tsne_and_tumor_extraction.pdf`
- `outputs/global_nonimmune_atlas/figures/global_nonimmune_tumor_cells_highlight_blue.pdf`
- `outputs/global_nonimmune_atlas/tables/global_nonimmune_counts.tsv`
- `outputs/global_nonimmune_atlas/tables/global_cluster_marker_summary.tsv`
- `outputs/global_nonimmune_atlas/tables/selected_global_tumor_clusters.tsv`
- `outputs/global_nonimmune_atlas/tables/all_global_tumor_summary.tsv`

What this step does:

- merges all processed matrices into one Seurat object
- computes a global t-SNE
- interprets clusters with broad marker programs
- identifies broad tumor-like clusters in the global map

## 2. Export a clean global t-SNE reference

Script:

- `scripts/export_global_tsne_reference_pdf.R`

This is only a formatting step used to keep a clean reference PDF with the same visual logic as downstream KCN plots.

## 3. Re-cluster all globally selected tumor cells with paper-like parameters

Script:

- `scripts/run_all_global_tumor_cells_exact_paper_pipeline.R`

Main outputs:

- `outputs/tumor_reclustering/tables/tumor_all_global_exact_summary.tsv`
- `outputs/tumor_reclustering/tables/tumor_all_global_cluster_sizes.tsv`
- `outputs/tumor_reclustering/tables/tumor_all_global_cluster_marker_summary.tsv`
- `outputs/tumor_reclustering/tables/tumor_all_global_cluster_sample_summary.tsv`

Paper-guided parameters used here:

- variable features with `mean.var.plot`
- `mean.cutoff = c(0.0125, 6)`
- `dispersion.cutoff = c(0.5, Inf)`
- UMI-regressed scaling
- `RunPCA(npcs = 30, seed.use = 12345)`
- clustering on PCs `1:20`
- `resolution = 0.8`
- `RunTSNE(dims = 1:5, seed.use = 12345)`

## 4. Annotate the tumor map into broad intestinal and diffuse poles

Script:

- `scripts/annotate_broad_tumor_tsne_intestinal_diffuse.R`

Main outputs:

- `outputs/tumor_poles/figures/tumor_all_global_intestinal_diffuse_annotated.pdf`
- `outputs/tumor_poles/tables/intestinal_diffuse_cell_labels.tsv`
- `outputs/tumor_poles/tables/intestinal_diffuse_cluster_summary.tsv`

Marker sets used:

- intestinal: `CDH17`, `REG4`, `MUC13`
- diffuse: `IGFBP5`, `COL1A1`, `S100A4`, `TAGLN`, `EGR1`

Interpretation:

- the goal here is not to recover every fine tumor subtype from the paper
- the goal is to recover a robust broad malignant axis useful for interpreting KCN genes

## 5. Project KCNQ1

Script:

- `scripts/plot_kcnq1_global_and_tumor.R`

Main outputs:

- `outputs/kcnq1/figures/KCNQ1_featureplot_global_and_tumor.pdf`
- `outputs/kcnq1/figures/KCNQ1_gastric_expression_tsne.pdf`
- `outputs/kcnq1/tables/KCNQ1_metrics_global_vs_tumor.tsv`

## 6. Project KCNE2 and KCNE3

Scripts:

- `scripts/plot_kcne2_kcne3_global_and_tumor.R`
- `scripts/run_kcne2_kcne3_spearman.R`
- `scripts/run_marker_spearman_check.R`

Main outputs:

- `outputs/kcne2_kcne3/figures/KCNE2_featureplot_global_and_tumor.pdf`
- `outputs/kcne2_kcne3/figures/KCNE3_featureplot_global_and_tumor.pdf`
- `outputs/kcne2_kcne3/tables/KCNE2_metrics_global_vs_tumor.tsv`
- `outputs/kcne2_kcne3/tables/KCNE3_metrics_global_vs_tumor.tsv`
- `outputs/kcne2_kcne3/tables/KCNE2_KCNE3_metrics_global_vs_tumor_combined.tsv`

## Main numbers

From the current outputs:

- global atlas: `13,022` cells
- broad tumor compartment: `1,408` cells
- re-clustered tumor subclusters: `12`

## Strengths of this branch

- easy to understand
- easy to rerun
- close to the published Seurat tumor-analysis logic
- enough to answer the main biological question at a broad level

## Limitations of this branch

- broad tumor selection is still marker-guided
- malignant calling is not fully reconstructed exactly as in the paper
- the final tumor structure is simplified into intestinal versus diffuse rather than every published substate

For a stricter robustness check, see `../pipeline_markers_plus_cnv/README.md`.
