suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Simple export focused on KCNQ1 in the final gastric pipeline.
project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
pipeline_dir <- file.path(project_dir, "pipeline_markers_only")

global_rds <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas", "objects", "global_nonimmune_seurat.rds")
tumor_rds <- file.path(
  pipeline_dir,
  "outputs",
  "tumor_reclustering",
  "objects",
  "tumor_all_global_exact_paper_seurat.rds"
)

output_dir <- file.path(pipeline_dir, "outputs", "kcnq1")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

global_obj <- readRDS(global_rds)
tumor_obj <- readRDS(tumor_rds)

if (!("KCNQ1" %in% rownames(global_obj)) || !("KCNQ1" %in% rownames(tumor_obj))) {
  stop("KCNQ1 is not present in one of the Seurat objects.")
}

compute_metrics <- function(obj, label) {
  raw_counts <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "counts")["KCNQ1", ])
  norm_values <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")["KCNQ1", ])
  positive_raw <- raw_counts[raw_counts > 0]

  data.frame(
    dataset = label,
    n_cells = length(raw_counts),
    n_positive_cells = length(positive_raw),
    pct_positive_cells = 100 * mean(raw_counts > 0),
    total_raw_transcripts = sum(raw_counts),
    mean_raw_counts_per_cell = mean(raw_counts),
    median_raw_counts_per_cell = median(raw_counts),
    mean_raw_counts_positive_cells = if (length(positive_raw) > 0) mean(positive_raw) else 0,
    median_raw_counts_positive_cells = if (length(positive_raw) > 0) median(positive_raw) else 0,
    max_raw_count = max(raw_counts),
    mean_log_normalized_expression = mean(norm_values),
    stringsAsFactors = FALSE
  )
}

metrics_df <- rbind(
  compute_metrics(global_obj, "Global non-immune t-SNE"),
  compute_metrics(tumor_obj, "Tumor-only t-SNE")
)

fwrite(metrics_df, file.path(table_dir, "KCNQ1_metrics_global_vs_tumor.tsv"), sep = "\t")
metrics_lines <- c(
  "KCNQ1 metrics in gastric single-cell analysis",
  "",
  "Global non-immune t-SNE",
  paste0("Cells: ", format(metrics_df$n_cells[1], big.mark = ",")),
  paste0("KCNQ1-positive cells: ", format(metrics_df$n_positive_cells[1], big.mark = ",")),
  paste0("% positive cells: ", sprintf("%.2f%%", metrics_df$pct_positive_cells[1])),
  paste0("Total raw transcripts: ", format(metrics_df$total_raw_transcripts[1], big.mark = ",")),
  paste0("Mean raw counts per cell: ", sprintf("%.3f", metrics_df$mean_raw_counts_per_cell[1])),
  paste0("Median raw counts per cell: ", sprintf("%.3f", metrics_df$median_raw_counts_per_cell[1])),
  paste0("Mean raw counts in positive cells: ", sprintf("%.3f", metrics_df$mean_raw_counts_positive_cells[1])),
  paste0("Median raw counts in positive cells: ", sprintf("%.3f", metrics_df$median_raw_counts_positive_cells[1])),
  paste0("Max raw count: ", format(metrics_df$max_raw_count[1], big.mark = ",")),
  paste0("Mean log-normalized expression: ", sprintf("%.3f", metrics_df$mean_log_normalized_expression[1])),
  "",
  "Tumor-only t-SNE",
  paste0("Cells: ", format(metrics_df$n_cells[2], big.mark = ",")),
  paste0("KCNQ1-positive cells: ", format(metrics_df$n_positive_cells[2], big.mark = ",")),
  paste0("% positive cells: ", sprintf("%.2f%%", metrics_df$pct_positive_cells[2])),
  paste0("Total raw transcripts: ", format(metrics_df$total_raw_transcripts[2], big.mark = ",")),
  paste0("Mean raw counts per cell: ", sprintf("%.3f", metrics_df$mean_raw_counts_per_cell[2])),
  paste0("Median raw counts per cell: ", sprintf("%.3f", metrics_df$median_raw_counts_per_cell[2])),
  paste0("Mean raw counts in positive cells: ", sprintf("%.3f", metrics_df$mean_raw_counts_positive_cells[2])),
  paste0("Median raw counts in positive cells: ", sprintf("%.3f", metrics_df$median_raw_counts_positive_cells[2])),
  paste0("Max raw count: ", format(metrics_df$max_raw_count[2], big.mark = ",")),
  paste0("Mean log-normalized expression: ", sprintf("%.3f", metrics_df$mean_log_normalized_expression[2]))
)
writeLines(metrics_lines, file.path(table_dir, "KCNQ1_metrics_global_vs_tumor.txt"))

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    aspect.ratio = 1
  )

global_plot <- FeaturePlot(
  global_obj,
  features = "KCNQ1",
  reduction = "tsne",
  cols = c("grey94", "#03045e"),
  pt.size = 0.45
) +
  labs(
    title = "KCNQ1 in the global non-immune t-SNE",
    subtitle = NULL
  ) +
  base_theme

tumor_plot <- FeaturePlot(
  tumor_obj,
  features = "KCNQ1",
  reduction = "tsne",
  cols = c("grey94", "#03045e"),
  pt.size = 0.8
) +
  labs(
    title = "KCNQ1 in the tumor-only t-SNE",
    subtitle = NULL
  ) +
  base_theme

combined_plot <- global_plot + tumor_plot

ggsave(
  file.path(figure_dir, "KCNQ1_featureplot_global_and_tumor.pdf"),
  combined_plot,
  width = 14,
  height = 7,
  device = cairo_pdf
)

message("Done. KCNQ1 feature plots and metrics exported.")
