suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Export a global t-SNE reference PDF using the same page size as the KCNE/KCNQ1 plots.
project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
pipeline_dir <- file.path(project_dir, "pipeline_markers_only")

global_rds <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas", "objects", "global_nonimmune_seurat.rds")
tumor_cells_tsv <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas", "tables", "all_global_tumor_cells.tsv")
output_pdf <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas", "figures", "global_nonimmune_tsne_reference_14x7.pdf")

global_obj <- readRDS(global_rds)
tumor_cells <- fread(tumor_cells_tsv)$cell_id

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

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
  cells.highlight = tumor_cells,
  cols = "grey85",
  cols.highlight = "#2b6cb0",
  pt.size = 0.45
) +
  labs(
    title = "Global non-immune t-SNE",
    subtitle = "Tumor-annotated cells highlighted",
    color = NULL
  ) +
  base_theme

ggsave(
  output_pdf,
  lineage_plot + tumor_highlight_plot,
  width = 14,
  height = 7,
  device = cairo_pdf
)

message("Done. Global t-SNE reference PDF exported.")
