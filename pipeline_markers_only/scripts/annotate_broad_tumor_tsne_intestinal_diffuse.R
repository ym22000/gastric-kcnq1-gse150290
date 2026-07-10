suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Step 3 of the final pipeline:
# annotate the broad tumor t-SNE with two simple poles.
# Intestinal pole = CDH17, REG4, MUC13
# Diffuse pole = IGFBP5, COL1A1, S100A4, TAGLN, EGR1
project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
pipeline_dir <- file.path(project_dir, "pipeline_markers_only")

input_rds <- file.path(
  pipeline_dir,
  "outputs",
  "tumor_reclustering",
  "objects",
  "tumor_all_global_exact_paper_seurat.rds"
)
output_dir <- file.path(pipeline_dir, "outputs", "tumor_poles")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)

intestinal_genes <- intersect(c("CDH17", "REG4", "MUC13"), rownames(obj))
diffuse_genes <- intersect(c("IGFBP5", "COL1A1", "S100A4", "TAGLN", "EGR1"), rownames(obj))

if (length(intestinal_genes) == 0 || length(diffuse_genes) == 0) {
  stop("Missing genes for one of the two pole signatures.")
}

obj <- AddModuleScore(
  obj,
  features = list(intestinal_genes),
  name = "IntestinalPole",
  assay = "RNA"
)
obj <- AddModuleScore(
  obj,
  features = list(diffuse_genes),
  name = "DiffusePole",
  assay = "RNA"
)

obj$IntestinalPoleScore <- obj$IntestinalPole1
obj$DiffusePoleScore <- obj$DiffusePole1
obj$PoleLabel <- ifelse(
  obj$IntestinalPoleScore >= obj$DiffusePoleScore,
  "Intestinal cells",
  "Diffuse cells"
)
obj$PoleLabel <- factor(obj$PoleLabel, levels = c("Intestinal cells", "Diffuse cells"))

cluster_summary <- obj@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    intestinal_score_mean = mean(IntestinalPoleScore),
    diffuse_score_mean = mean(DiffusePoleScore),
    intestinal_fraction = mean(PoleLabel == "Intestinal cells"),
    diffuse_fraction = mean(PoleLabel == "Diffuse cells"),
    assigned_pole = ifelse(mean(PoleLabel == "Intestinal cells") >= 0.5, "Intestinal cells", "Diffuse cells"),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cell_summary <- obj@meta.data %>%
  mutate(
    cell_id = rownames(obj@meta.data),
    cluster = as.character(seurat_clusters)
  ) %>%
  select(cell_id, cluster, sample_title, tissue_group, IntestinalPoleScore, DiffusePoleScore, PoleLabel)

fwrite(cluster_summary, file.path(table_dir, "intestinal_diffuse_cluster_summary.tsv"), sep = "\t")
fwrite(cell_summary, file.path(table_dir, "intestinal_diffuse_cell_labels.tsv"), sep = "\t")

pole_cols <- c(
  "Intestinal cells" = "#d7301f",
  "Diffuse cells" = "#7b3294"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    aspect.ratio = 1
  )

state_plot <- DimPlot(
  obj,
  reduction = "tsne",
  group.by = "PoleLabel",
  cols = pole_cols,
  pt.size = 0.8
) +
  labs(
    title = "Tumor cells",
    subtitle = "Intestinal versus diffuse poles",
    color = NULL
  ) +
  base_theme

marker_genes <- c("CDH17", "REG4", "MUC13", "IGFBP5", "COL1A1", "S100A4", "TAGLN", "EGR1")
marker_genes <- intersect(marker_genes, rownames(obj))

marker_plots <- lapply(marker_genes, function(g) {
  FeaturePlot(
    obj,
    features = g,
    reduction = "tsne",
    cols = c("grey94", "#c62828"),
    pt.size = 0.8
  ) +
    labs(title = g, subtitle = NULL) +
    base_theme +
    theme(legend.position = "right")
})

all_plots <- c(list(state_plot), marker_plots)
combined <- wrap_plots(all_plots, ncol = 3)

ggsave(
  file.path(figure_dir, "tumor_all_global_intestinal_diffuse_annotated.pdf"),
  combined,
  width = 14,
  height = 14,
  device = cairo_pdf
)

message("Done. Intestinal versus diffuse annotated t-SNE exported.")
