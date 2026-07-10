suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
marker_pipeline_dir <- file.path(project_dir, "pipeline_markers_only")
cnv_pipeline_dir <- file.path(project_dir, "pipeline_markers_plus_cnv")

global_rds <- file.path(
  marker_pipeline_dir,
  "outputs",
  "global_nonimmune_atlas",
  "objects",
  "global_nonimmune_seurat.rds"
)
tumor_rds <- file.path(
  marker_pipeline_dir,
  "outputs",
  "tumor_reclustering",
  "objects",
  "tumor_all_global_exact_paper_seurat.rds"
)
full_gene_rds <- file.path(
  cnv_pipeline_dir,
  "outputs",
  "infercnv_validation",
  "objects",
  "selected_cells_full_gene_space_seurat.rds"
)
infercnv_rds <- file.path(
  cnv_pipeline_dir,
  "outputs",
  "infercnv_validation",
  "objects",
  "infercnv_result_object.rds"
)
refined_rds <- file.path(
  cnv_pipeline_dir,
  "outputs",
  "tumor_refined_markers_cnv",
  "objects",
  "tumor_refined_markers_cnv_seurat.rds"
)

step1_script <- file.path(marker_pipeline_dir, "scripts", "run_global_nonimmune_tumor_extraction.R")
step2_script <- file.path(marker_pipeline_dir, "scripts", "run_all_global_tumor_cells_exact_paper_pipeline.R")
infercnv_script <- file.path(cnv_pipeline_dir, "scripts", "run_infercnv_on_existing_global_tsne.R")
refine_script <- file.path(cnv_pipeline_dir, "scripts", "refine_tumor_compartment_markers_plus_cnv.R")

run_script <- function(path) {
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
}

if (!file.exists(global_rds)) {
  message("Global Seurat object missing. Rebuilding step 1...")
  run_script(step1_script)
}

if (!file.exists(tumor_rds)) {
  message("Tumor Seurat object missing. Rebuilding step 2...")
  run_script(step2_script)
}

if (!file.exists(full_gene_rds) || !file.exists(infercnv_rds)) {
  message("inferCNV objects missing. Rebuilding inferCNV validation step...")
  run_script(infercnv_script)
}

if (!file.exists(refined_rds)) {
  message("Refined tumor object missing. Rebuilding refinement step...")
  run_script(refine_script)
}

output_dir <- file.path(cnv_pipeline_dir, "outputs", "refined_tumor_marker_sets_and_kcn")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(refined_rds)

if (!("PoleLabel" %in% colnames(obj@meta.data))) {
  intestinal_genes <- intersect(c("CDH17", "REG4", "MUC13"), rownames(obj))
  diffuse_genes <- intersect(c("IGFBP5", "COL1A1", "S100A4", "TAGLN", "EGR1"), rownames(obj))
  obj <- AddModuleScore(obj, features = list(intestinal_genes), name = "IntestinalPole", assay = "RNA")
  obj <- AddModuleScore(obj, features = list(diffuse_genes), name = "DiffusePole", assay = "RNA")
  obj$IntestinalPoleScore <- obj$IntestinalPole1
  obj$DiffusePoleScore <- obj$DiffusePole1
  obj$PoleLabel <- ifelse(
    obj$IntestinalPoleScore >= obj$DiffusePoleScore,
    "Intestinal cells",
    "Diffuse cells"
  )
}
obj$PoleLabel <- factor(obj$PoleLabel, levels = c("Intestinal cells", "Diffuse cells"))

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

make_reference_plot <- function() {
  DimPlot(
    obj,
    reduction = "tsne",
    group.by = "PoleLabel",
    cols = c("Intestinal cells" = "#d7301f", "Diffuse cells" = "#7b3294"),
    pt.size = 0.8
  ) +
    labs(
      title = "Refined tumor-only t-SNE",
      subtitle = NULL,
      color = NULL
    ) +
    base_theme
}

make_feature_plot <- function(gene, color = "#c62828") {
  FeaturePlot(
    obj,
    features = gene,
    reduction = "tsne",
    cols = c("grey94", color),
    pt.size = 0.8
  ) +
    labs(title = gene, subtitle = NULL) +
    base_theme
}

make_score_plot <- function(feature_name, title_text, color = "#c62828") {
  FeaturePlot(
    obj,
    features = feature_name,
    reduction = "tsne",
    cols = c("grey94", color),
    pt.size = 0.8
  ) +
    labs(title = title_text, subtitle = NULL) +
    base_theme
}

pdf(
  file.path(figure_dir, "refined_tumor_only_tsne_cnv_intestinal_diffuse_kcn_multiplot.pdf"),
  width = 10,
  height = 8,
  onefile = TRUE
)

print(wrap_plots(
  make_reference_plot(),
  make_feature_plot("CD44"),
  make_feature_plot("CDK6"),
  make_feature_plot("GATA4"),
  ncol = 2
))

print(wrap_plots(
  make_reference_plot(),
  make_feature_plot("GATA6"),
  make_feature_plot("KLF5"),
  make_feature_plot("KRAS"),
  ncol = 2
))

print(wrap_plots(
  make_reference_plot(),
  make_feature_plot("CDH17"),
  make_feature_plot("REG4"),
  make_feature_plot("MUC13"),
  ncol = 2
))

print(wrap_plots(
  make_reference_plot(),
  make_feature_plot("IGFBP5"),
  make_feature_plot("COL1A1"),
  make_feature_plot("S100A4"),
  ncol = 2
))

print(wrap_plots(
  make_reference_plot(),
  make_feature_plot("TAGLN"),
  make_feature_plot("EGR1"),
  make_score_plot("DiffusePoleScore", "Diffuse score"),
  ncol = 2
))

print(wrap_plots(
  make_reference_plot(),
  make_feature_plot("KCNQ1", color = "#03045e"),
  make_feature_plot("KCNE2", color = "#03045e"),
  make_feature_plot("KCNE3", color = "#03045e"),
  ncol = 2
))

dev.off()

message("Done. Refined tumor-only marker and KCN multiplot PDF exported.")
