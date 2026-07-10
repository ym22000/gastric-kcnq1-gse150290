suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Export the global non-immune t-SNE together with tumor markers and
# CNV-like genes discussed in the gastric paper.
project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
marker_pipeline_dir <- file.path(project_dir, "pipeline_markers_only")
cnv_pipeline_dir <- file.path(project_dir, "pipeline_markers_plus_cnv")

global_rds <- file.path(marker_pipeline_dir, "outputs", "global_nonimmune_atlas", "objects", "global_nonimmune_seurat.rds")
tumor_cells_tsv <- file.path(marker_pipeline_dir, "outputs", "global_nonimmune_atlas", "tables", "all_global_tumor_cells.tsv")

output_dir <- file.path(cnv_pipeline_dir, "outputs", "global_cnv_screen")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

global_obj <- readRDS(global_rds)
tumor_cells <- fread(tumor_cells_tsv)$cell_id

tumor_marker_genes <- c("EPCAM", "CDH17", "COL3A1", "PDGFRB")
cnv_like_genes <- c("CD44", "CDK6", "GATA4", "GATA6", "KLF5", "KRAS")
all_genes <- c(tumor_marker_genes, cnv_like_genes)

missing_genes <- setdiff(all_genes, rownames(global_obj))
if (length(missing_genes) > 0) {
  stop(
    paste0(
      "The following genes are missing from the global object: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

global_obj$tumor_screen_group <- ifelse(
  colnames(global_obj) %in% tumor_cells,
  "Extracted tumor cells",
  "Other non-immune cells"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

make_feature_plot <- function(gene) {
  FeaturePlot(
    global_obj,
    features = gene,
    reduction = "tsne",
    cols = c("grey94", "#b2182b"),
    pt.size = 0.45
  ) +
    labs(
      title = gene,
      subtitle = NULL
    ) +
    base_theme
}

lineage_plot <- DimPlot(
  global_obj,
  reduction = "tsne",
  group.by = "predicted_lineage",
  label = FALSE,
  pt.size = 0.45
) +
  labs(
    title = "Global non-immune t-SNE",
    subtitle = "Marker-based broad lineages",
    color = "Lineage"
  ) +
  base_theme

tumor_highlight_plot <- DimPlot(
  global_obj,
  reduction = "tsne",
  group.by = "tumor_screen_group",
  cols = c("Other non-immune cells" = "grey85", "Extracted tumor cells" = "#2b6cb0"),
  pt.size = 0.45
) +
  labs(
    title = "Global non-immune t-SNE",
    subtitle = "Current extracted tumor compartment",
    color = NULL
  ) +
  base_theme

tumor_marker_plots <- lapply(tumor_marker_genes, make_feature_plot)
cnv_like_plots <- lapply(cnv_like_genes, make_feature_plot)

pdf(
  file.path(figure_dir, "global_tsne_tumor_markers_and_cnv_like_genes.pdf"),
  width = 14,
  height = 7,
  onefile = TRUE
)

print(lineage_plot + tumor_highlight_plot)
print(wrap_plots(tumor_marker_plots, ncol = 2))
print(wrap_plots(cnv_like_plots[1:4], ncol = 2))
print(wrap_plots(cnv_like_plots[5:6], ncol = 2))

dev.off()

compute_gene_metrics <- function(obj, gene, group_var) {
  raw_counts <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "counts")[gene, ])
  norm_values <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")[gene, ])
  meta <- obj@meta.data

  dt <- data.table(
    gene = gene,
    group = as.character(meta[[group_var]]),
    raw_count = raw_counts,
    norm_value = norm_values
  )

  dt[
    ,
    .(
      n_cells = .N,
      n_positive_cells = sum(raw_count > 0),
      pct_positive_cells = 100 * mean(raw_count > 0),
      total_raw_transcripts = sum(raw_count),
      mean_raw_counts_per_cell = mean(raw_count),
      mean_log_normalized_expression = mean(norm_value),
      max_raw_count = max(raw_count)
    ),
    by = .(gene, group)
  ]
}

screen_table <- rbindlist(lapply(all_genes, compute_gene_metrics, obj = global_obj, group_var = "tumor_screen_group"))
lineage_table <- rbindlist(lapply(all_genes, compute_gene_metrics, obj = global_obj, group_var = "predicted_lineage"))

fwrite(screen_table, file.path(table_dir, "tumor_screen_gene_metrics.tsv"), sep = "\t")
fwrite(lineage_table, file.path(table_dir, "predicted_lineage_gene_metrics.tsv"), sep = "\t")

message("Done. Global tumor/CNV screening plots and tables exported.")
