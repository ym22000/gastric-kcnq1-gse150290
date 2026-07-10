# Pipeline Markers Plus CNV

This branch is the robustness extension of the `GSE150290` gastric analysis.

It starts from the same dataset and the same broad tumor-oriented logic as `pipeline_markers_only`, but it adds CNV-informed filtering to ask whether the interpretation of `KCNQ1`, `KCNE2`, and `KCNE3` remains stable in a stricter malignant compartment.

## Input

Shared upstream context:

- `../cache/`
- `../raw_geo/`
- `../pipeline_markers_only/outputs/`

Main reference documents:

- `../bibliography/article_gastrique.pdf`
- `../bibliography/supplementary_data_gastrique.pdf`

## Script order

### 1. Visual CNV-oriented screen on the global map

Script:

- `scripts/plot_global_tumor_and_cnv_markers.R`

Main outputs:

- `outputs/global_cnv_screen/figures/global_tsne_tumor_markers_and_cnv_like_genes.pdf`
- `outputs/global_cnv_screen/tables/tumor_screen_gene_metrics.tsv`
- `outputs/global_cnv_screen/tables/predicted_lineage_gene_metrics.tsv`

Purpose:

- project epithelial tumor markers
- project CNV-like or malignancy-associated genes discussed in the paper
- visually check whether the broad tumor compartment already captures the strongest malignant signal

### 2. Run inferCNV on the already selected cells

Script:

- `scripts/run_infercnv_on_existing_global_tsne.R`

Main outputs kept in the cleaned project:

- `outputs/infercnv_validation/figures/infercnv_validation_global_and_tumor_tsne.pdf`
- `outputs/infercnv_validation/tables/infercnv_cell_annotations.tsv`
- `outputs/infercnv_validation/tables/infercnv_cell_scores.tsv`
- `outputs/infercnv_validation/tables/infercnv_reference_cells.tsv`
- `outputs/infercnv_validation/tables/infercnv_reference_summary.tsv`
- `outputs/infercnv_validation/tables/infercnv_tumor_enrichment.tsv`
- `outputs/infercnv_validation/tables/hg19_gene_order_for_infercnv.tsv`

Note:

Very large inferCNV intermediate state files were removed from this cleaned repository because they are fully regenerable from the scripts and raw inputs.

### 3. Refine the tumor compartment with markers plus CNV support

Script:

- `scripts/refine_tumor_compartment_markers_plus_cnv.R`

Main outputs:

- `outputs/tumor_refined_markers_cnv/figures/refined_tumor_markers_plus_cnv_validation.pdf`
- `outputs/tumor_refined_markers_cnv/tables/current_vs_refined_tumor_comparison.tsv`
- `outputs/tumor_refined_markers_cnv/tables/refined_added_cells_by_lineage.tsv`
- `outputs/tumor_refined_markers_cnv/tables/refined_tumor_cells.tsv`
- `outputs/tumor_refined_markers_cnv/tables/refined_tumor_cluster_sizes.tsv`
- `outputs/tumor_refined_markers_cnv/tables/refined_tumor_pole_labels.tsv`
- `outputs/tumor_refined_markers_cnv/tables/refined_tumor_summary.tsv`

Purpose:

- keep cells already supported by the broad marker-based tumor call
- add cells with compatible epithelial identity and stronger CNV support
- avoid using CNV alone as the only malignancy criterion

### 4. Re-evaluate the KCN genes in the refined tumor compartment

Script:

- `scripts/plot_refined_kcnq1_kcne2_kcne3_and_spearman.R`

Main outputs:

- `outputs/refined_kcnq1_kcne2_kcne3/figures/refined_tumor_KCNQ1_KCNE2_KCNE3_multiplot.pdf`
- `outputs/refined_kcnq1_kcne2_kcne3/tables/refined_KCNQ1_KCNE2_KCNE3_spearman_metrics.tsv`

## Main numbers

From `outputs/tumor_refined_markers_cnv/tables/refined_tumor_summary.tsv`:

- current tumor compartment: `1,402` cells
- refined tumor compartment: `1,614` cells
- retained from current compartment: `1,083` cells
- added by markers plus CNV support: `531` cells
- CNV-high fraction rises from `54.1%` to `79.9%`

## Why this branch matters

This branch is useful because it asks a more rigorous question:

If the marker-based tumor compartment is imperfect, do the KCN conclusions still hold after adding a second orthogonal filter?

In practice, the answer is yes:

- `KCNQ1` remains intestinal-oriented
- `KCNE3` remains intestinal-oriented
- `KCNE2` remains weak

## What this branch is and is not

This branch is:

- a robustness analysis on the same dataset
- a way to strengthen the malignant interpretation
- a protection against over-interpreting marker-only selection

This branch is not:

- a strict reproduction of every malignancy call from the original paper
- a standalone proof that CNV alone defines gastric tumor cells

The correct interpretation is therefore:

Markers define the broad tumor context, CNV adds support, and the KCN conclusions are considered stronger only when both point in the same direction.
