suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
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

full_gene_rds <- file.path(
  cnv_pipeline_dir,
  "outputs",
  "infercnv_validation",
  "objects",
  "selected_cells_full_gene_space_seurat.rds"
)
cnv_tsv <- file.path(
  cnv_pipeline_dir,
  "outputs",
  "infercnv_validation",
  "tables",
  "infercnv_cell_scores.tsv"
)
current_tumor_rds <- file.path(
  marker_pipeline_dir,
  "outputs",
  "tumor_reclustering",
  "objects",
  "tumor_all_global_exact_paper_seurat.rds"
)
current_poles_tsv <- file.path(
  marker_pipeline_dir,
  "outputs",
  "tumor_poles",
  "tables",
  "intestinal_diffuse_cell_labels.tsv"
)

output_dir <- file.path(cnv_pipeline_dir, "outputs", "tumor_refined_markers_cnv")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

refined_rds <- file.path(object_dir, "tumor_refined_markers_cnv_seurat.rds")

full_obj <- readRDS(full_gene_rds)
current_tumor_obj <- readRDS(current_tumor_rds)
current_poles <- fread(current_poles_tsv)
cnv <- fread(cnv_tsv)

meta <- as.data.table(full_obj@meta.data)
meta[, cell_id := rownames(full_obj@meta.data)]
meta <- merge(meta, cnv[, .(cell_id, cnv_score, cnv_high, in_current_tumor_compartment)], by = "cell_id", all.x = TRUE)

score_cells <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    return(rep(0, ncol(obj)))
  }
  Matrix::colMeans(GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE])
}

meta[, epithelial_score2 := score_cells(full_obj, c("EPCAM", "KRT19", "KRT8", "KRT18", "KRT17", "MUC1", "TFF3", "TFF2", "MUC6", "CDH17"))]
meta[, stromal_score2 := score_cells(full_obj, c("COL1A1", "COL3A1", "TAGLN", "DCN", "PDGFRB", "PDGFRA", "ACTA2", "THY1"))]
meta[, immune_score2 := score_cells(full_obj, c("PTPRC", "LST1", "FCN1", "AIF1", "TYMP", "CTSS", "HLA-DRA", "LYZ", "S100A8", "S100A9"))]
meta[, tumor_program_score := score_cells(full_obj, c("EPCAM", "CDH17", "COL3A1", "PDGFRB"))]

current_tumor_meta <- meta[in_current_tumor_compartment == TRUE]
epi_q75 <- as.numeric(quantile(current_tumor_meta$epithelial_score2, 0.75, na.rm = TRUE))
tumor_program_med <- median(current_tumor_meta$tumor_program_score, na.rm = TRUE)

# Conservative refined rule:
# 1. Keep current tumor cells if they are CNV-high OR if they retain at least the
#    median tumor-program score of the current compartment.
# 2. Add new cells only if they are CNV-high, come from cancer lesions, belong to
#    epithelial-like global lineages, and have epithelial dominance over stromal
#    and immune programs.
meta[, refined_rule := "discarded"]
meta[in_current_tumor_compartment == TRUE & (cnv_high == TRUE | tumor_program_score >= tumor_program_med), refined_rule := "kept_from_current_tumor"]
meta[
  in_current_tumor_compartment == FALSE &
    cnv_high == TRUE &
    tissue_group == "Cancer" &
    predicted_lineage %in% c("GMC", "IM", "PMC", "enteroendocrine") &
    epithelial_score2 >= epi_q75 &
    epithelial_score2 > stromal_score2 &
    epithelial_score2 > immune_score2,
  refined_rule := "added_by_markers_plus_cnv"
]

refined_cells <- meta[refined_rule != "discarded", cell_id]
meta[, refined_tumor_compartment := cell_id %in% refined_cells]

summary_dt <- data.table(
  metric = c(
    "current_tumor_cells",
    "refined_tumor_cells",
    "kept_from_current_tumor",
    "added_by_markers_plus_cnv",
    "current_tumor_cnv_high_pct",
    "refined_tumor_cnv_high_pct",
    "epi_q75_current_tumor",
    "tumor_program_median_current_tumor"
  ),
  value = c(
    sum(meta$in_current_tumor_compartment, na.rm = TRUE),
    length(refined_cells),
    sum(meta$refined_rule == "kept_from_current_tumor"),
    sum(meta$refined_rule == "added_by_markers_plus_cnv"),
    100 * mean(meta$cnv_high[meta$in_current_tumor_compartment], na.rm = TRUE),
    100 * mean(meta$cnv_high[meta$refined_tumor_compartment], na.rm = TRUE),
    epi_q75,
    tumor_program_med
  )
)
fwrite(summary_dt, file.path(table_dir, "refined_tumor_summary.tsv"), sep = "\t")

rule_cells_dt <- meta[
  refined_rule != "discarded",
  .(
    cell_id,
    refined_rule,
    predicted_lineage,
    sample_title,
    tissue_group,
    cnv_score,
    cnv_high,
    epithelial_score2,
    stromal_score2,
    immune_score2,
    tumor_program_score
  )
]
fwrite(rule_cells_dt, file.path(table_dir, "refined_tumor_cells.tsv"), sep = "\t")

added_by_lineage <- meta[refined_rule == "added_by_markers_plus_cnv", .N, by = predicted_lineage][order(-N)]
fwrite(added_by_lineage, file.path(table_dir, "refined_added_cells_by_lineage.tsv"), sep = "\t")

counts_refined <- GetAssayData(full_obj, assay = "RNA", layer = "counts")[, refined_cells, drop = FALSE]
meta_refined <- as.data.frame(full_obj@meta.data[refined_cells, , drop = FALSE])
meta_refined$refined_rule <- meta$refined_rule[match(refined_cells, meta$cell_id)]
meta_refined$cnv_score <- meta$cnv_score[match(refined_cells, meta$cell_id)]
meta_refined$cnv_high <- meta$cnv_high[match(refined_cells, meta$cell_id)]

run_paper_like_seurat <- function(seu) {
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(
    seu,
    selection.method = "mean.var.plot",
    mean.cutoff = c(0.0125, 6),
    dispersion.cutoff = c(0.5, Inf),
    verbose = FALSE
  )
  seu <- ScaleData(
    seu,
    features = VariableFeatures(seu),
    vars.to.regress = "nCount_RNA",
    verbose = FALSE
  )
  seu <- RunPCA(
    seu,
    features = VariableFeatures(seu),
    npcs = 30,
    seed.use = 12345,
    verbose = FALSE
  )
  seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)
  seu <- RunTSNE(
    seu,
    dims = 1:5,
    seed.use = 12345,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  seu
}

refined_obj <- CreateSeuratObject(
  counts = counts_refined,
  meta.data = meta_refined,
  project = "GSE150290_refined_tumor_markers_cnv"
)
refined_obj <- run_paper_like_seurat(refined_obj)

intestinal_genes <- intersect(c("CDH17", "REG4", "MUC13"), rownames(refined_obj))
diffuse_genes <- intersect(c("IGFBP5", "COL1A1", "S100A4", "TAGLN", "EGR1"), rownames(refined_obj))

refined_obj <- AddModuleScore(refined_obj, features = list(intestinal_genes), name = "IntestinalPole", assay = "RNA")
refined_obj <- AddModuleScore(refined_obj, features = list(diffuse_genes), name = "DiffusePole", assay = "RNA")
refined_obj$IntestinalPoleScore <- refined_obj$IntestinalPole1
refined_obj$DiffusePoleScore <- refined_obj$DiffusePole1
refined_obj$PoleLabel <- ifelse(
  refined_obj$IntestinalPoleScore >= refined_obj$DiffusePoleScore,
  "Intestinal cells",
  "Diffuse cells"
)
refined_obj$PoleLabel <- factor(refined_obj$PoleLabel, levels = c("Intestinal cells", "Diffuse cells"))

cluster_sizes <- as.data.frame(table(refined_obj$seurat_clusters), stringsAsFactors = FALSE)
colnames(cluster_sizes) <- c("cluster", "n_cells")
fwrite(cluster_sizes, file.path(table_dir, "refined_tumor_cluster_sizes.tsv"), sep = "\t")

refined_cell_labels <- refined_obj@meta.data %>%
  mutate(cell_id = rownames(refined_obj@meta.data)) %>%
  select(cell_id, refined_rule, sample_title, tissue_group, cnv_score, cnv_high, IntestinalPoleScore, DiffusePoleScore, PoleLabel, seurat_clusters)
fwrite(refined_cell_labels, file.path(table_dir, "refined_tumor_pole_labels.tsv"), sep = "\t")

comparison_dt <- data.table(
  version = c("Current tumor map", "Refined markers+CNV map"),
  n_cells = c(ncol(current_tumor_obj), ncol(refined_obj)),
  intestinal_cells = c(
    sum(current_poles$PoleLabel == "Intestinal cells"),
    sum(refined_obj$PoleLabel == "Intestinal cells")
  ),
  diffuse_cells = c(
    sum(current_poles$PoleLabel == "Diffuse cells"),
    sum(refined_obj$PoleLabel == "Diffuse cells")
  ),
  mean_cnv_score = c(
    mean(cnv$cnv_score[match(colnames(current_tumor_obj), cnv$cell_id)], na.rm = TRUE),
    mean(refined_obj$cnv_score, na.rm = TRUE)
  ),
  pct_cnv_high = c(
    100 * mean(cnv$cnv_high[match(colnames(current_tumor_obj), cnv$cell_id)], na.rm = TRUE),
    100 * mean(refined_obj$cnv_high, na.rm = TRUE)
  )
)
fwrite(comparison_dt, file.path(table_dir, "current_vs_refined_tumor_comparison.tsv"), sep = "\t")

saveRDS(refined_obj, refined_rds)

# Plotting
base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

full_obj$refined_group <- ifelse(
  colnames(full_obj) %in% refined_cells,
  "Refined tumor cells",
  "Other selected cells"
)
full_obj$current_group <- ifelse(
  colnames(full_obj) %in% current_poles$cell_id,
  "Current tumor cells",
  "Other selected cells"
)
full_obj$cnv_high_label <- ifelse(
  cnv$cnv_high[match(colnames(full_obj), cnv$cell_id)],
  "CNV-high",
  "CNV-low"
)
full_obj$cnv_score <- cnv$cnv_score[match(colnames(full_obj), cnv$cell_id)]

current_global_plot <- DimPlot(
  full_obj,
  reduction = "tsne",
  group.by = "current_group",
  cols = c("Other selected cells" = "grey85", "Current tumor cells" = "#2b6cb0"),
  pt.size = 0.45
) +
  labs(title = "Global gastric t-SNE", subtitle = "Current tumor compartment", color = NULL) +
  base_theme

refined_global_plot <- DimPlot(
  full_obj,
  reduction = "tsne",
  group.by = "refined_group",
  cols = c("Other selected cells" = "grey85", "Refined tumor cells" = "#d73027"),
  pt.size = 0.45
) +
  labs(title = "Global gastric t-SNE", subtitle = "Refined markers + CNV tumor compartment", color = NULL) +
  base_theme

global_cnv_plot <- DimPlot(
  full_obj,
  reduction = "tsne",
  group.by = "cnv_high_label",
  cols = c("CNV-low" = "grey85", "CNV-high" = "#7b2cbf"),
  pt.size = 0.45
) +
  labs(title = "Global gastric t-SNE", subtitle = "CNV-high cells", color = NULL) +
  base_theme

current_tumor_obj$PoleLabel <- current_poles$PoleLabel[match(colnames(current_tumor_obj), current_poles$cell_id)]
current_tumor_obj$PoleLabel <- factor(current_tumor_obj$PoleLabel, levels = c("Intestinal cells", "Diffuse cells"))
current_tumor_obj$cnv_score <- cnv$cnv_score[match(colnames(current_tumor_obj), cnv$cell_id)]

pole_cols <- c("Intestinal cells" = "#d7301f", "Diffuse cells" = "#7b3294")

current_pole_plot <- DimPlot(
  current_tumor_obj,
  reduction = "tsne",
  group.by = "PoleLabel",
  cols = pole_cols,
  pt.size = 0.8
) +
  labs(title = "Tumor-only t-SNE", subtitle = "Current map", color = NULL) +
  base_theme

refined_pole_plot <- DimPlot(
  refined_obj,
  reduction = "tsne",
  group.by = "PoleLabel",
  cols = pole_cols,
  pt.size = 0.8
) +
  labs(title = "Tumor-only t-SNE", subtitle = "Refined markers + CNV map", color = NULL) +
  base_theme

make_refined_feature <- function(gene, color = "#c62828") {
  FeaturePlot(
    refined_obj,
    features = gene,
    reduction = "tsne",
    cols = c("grey94", color),
    pt.size = 0.8
  ) +
    labs(title = gene, subtitle = NULL) +
    base_theme +
    theme(legend.position = "right")
}

intestinal_marker_plots <- lapply(intersect(c("CDH17", "REG4", "MUC13"), rownames(refined_obj)), make_refined_feature)
diffuse_marker_plots <- lapply(intersect(c("IGFBP5", "COL1A1", "S100A4", "TAGLN", "EGR1"), rownames(refined_obj)), make_refined_feature)
cnv_marker_plots <- lapply(intersect(c("CD44", "CDK6", "GATA4", "GATA6", "KLF5", "KRAS"), rownames(refined_obj)), make_refined_feature, color = "#7b2cbf")

pdf(
  file.path(figure_dir, "refined_tumor_markers_plus_cnv_validation.pdf"),
  width = 14,
  height = 7,
  onefile = TRUE
)

print(current_global_plot + refined_global_plot + global_cnv_plot + plot_layout(ncol = 3))
print(current_pole_plot + refined_pole_plot)
print(wrap_plots(c(list(refined_pole_plot), intestinal_marker_plots), ncol = 2))
print(wrap_plots(c(list(refined_pole_plot), diffuse_marker_plots), ncol = 3))
print(wrap_plots(c(list(refined_pole_plot), cnv_marker_plots), ncol = 3))

dev.off()

message("Done. Refined tumor compartment exported.")
