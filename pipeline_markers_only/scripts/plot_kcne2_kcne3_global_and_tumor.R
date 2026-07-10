suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Export KCNE2 and KCNE3 on the final gastric global and tumor t-SNE objects.
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

output_dir <- file.path(pipeline_dir, "outputs", "kcne2_kcne3")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

global_obj <- readRDS(global_rds)
tumor_obj <- readRDS(tumor_rds)

target_genes <- c("KCNE2", "KCNE3")
missing_global <- setdiff(target_genes, rownames(global_obj))
missing_tumor <- setdiff(target_genes, rownames(tumor_obj))

if (length(missing_global) > 0 || length(missing_tumor) > 0) {
  stop(
    paste(
      "Missing genes.",
      if (length(missing_global) > 0) {
        paste0("Global object: ", paste(missing_global, collapse = ", "))
      } else {
        NULL
      },
      if (length(missing_tumor) > 0) {
        paste0("Tumor object: ", paste(missing_tumor, collapse = ", "))
      } else {
        NULL
      }
    )
  )
}

compute_metrics <- function(obj, gene, label) {
  raw_counts <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "counts")[gene, ])
  norm_values <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")[gene, ])
  positive_raw <- raw_counts[raw_counts > 0]

  data.frame(
    gene = gene,
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

make_metrics_text <- function(metrics_df, gene) {
  global_row <- metrics_df[metrics_df$dataset == "Global non-immune t-SNE", ]
  tumor_row <- metrics_df[metrics_df$dataset == "Tumor-only t-SNE", ]

  c(
    paste0(gene, " metrics in gastric single-cell analysis"),
    "",
    "Global non-immune t-SNE",
    paste0("Cells: ", format(global_row$n_cells, big.mark = ",")),
    paste0(gene, "-positive cells: ", format(global_row$n_positive_cells, big.mark = ",")),
    paste0("% positive cells: ", sprintf("%.2f%%", global_row$pct_positive_cells)),
    paste0("Total raw transcripts: ", format(global_row$total_raw_transcripts, big.mark = ",")),
    paste0("Mean raw counts per cell: ", sprintf("%.3f", global_row$mean_raw_counts_per_cell)),
    paste0("Median raw counts per cell: ", sprintf("%.3f", global_row$median_raw_counts_per_cell)),
    paste0("Mean raw counts in positive cells: ", sprintf("%.3f", global_row$mean_raw_counts_positive_cells)),
    paste0("Median raw counts in positive cells: ", sprintf("%.3f", global_row$median_raw_counts_positive_cells)),
    paste0("Max raw count: ", format(global_row$max_raw_count, big.mark = ",")),
    paste0("Mean log-normalized expression: ", sprintf("%.3f", global_row$mean_log_normalized_expression)),
    "",
    "Tumor-only t-SNE",
    paste0("Cells: ", format(tumor_row$n_cells, big.mark = ",")),
    paste0(gene, "-positive cells: ", format(tumor_row$n_positive_cells, big.mark = ",")),
    paste0("% positive cells: ", sprintf("%.2f%%", tumor_row$pct_positive_cells)),
    paste0("Total raw transcripts: ", format(tumor_row$total_raw_transcripts, big.mark = ",")),
    paste0("Mean raw counts per cell: ", sprintf("%.3f", tumor_row$mean_raw_counts_per_cell)),
    paste0("Median raw counts per cell: ", sprintf("%.3f", tumor_row$median_raw_counts_per_cell)),
    paste0("Mean raw counts in positive cells: ", sprintf("%.3f", tumor_row$mean_raw_counts_positive_cells)),
    paste0("Median raw counts in positive cells: ", sprintf("%.3f", tumor_row$median_raw_counts_positive_cells)),
    paste0("Max raw count: ", format(tumor_row$max_raw_count, big.mark = ",")),
    paste0("Mean log-normalized expression: ", sprintf("%.3f", tumor_row$mean_log_normalized_expression))
  )
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    aspect.ratio = 1
  )

all_metrics <- lapply(target_genes, function(gene) {
  metrics_df <- rbind(
    compute_metrics(global_obj, gene, "Global non-immune t-SNE"),
    compute_metrics(tumor_obj, gene, "Tumor-only t-SNE")
  )

  fwrite(
    metrics_df,
    file.path(table_dir, paste0(gene, "_metrics_global_vs_tumor.tsv")),
    sep = "\t"
  )
  writeLines(
    make_metrics_text(metrics_df, gene),
    file.path(table_dir, paste0(gene, "_metrics_global_vs_tumor.txt"))
  )

  global_plot <- FeaturePlot(
    global_obj,
    features = gene,
    reduction = "tsne",
    cols = c("grey94", "#03045e"),
    pt.size = 0.45
  ) +
    labs(title = paste0(gene, " in the global non-immune t-SNE"), subtitle = NULL) +
    base_theme

  tumor_plot <- FeaturePlot(
    tumor_obj,
    features = gene,
    reduction = "tsne",
    cols = c("grey94", "#03045e"),
    pt.size = 0.8
  ) +
    labs(title = paste0(gene, " in the tumor-only t-SNE"), subtitle = NULL) +
    base_theme

  ggsave(
    file.path(figure_dir, paste0(gene, "_featureplot_global_and_tumor.pdf")),
    global_plot + tumor_plot,
    width = 14,
    height = 7,
    device = cairo_pdf
  )

  metrics_df
})

combined_metrics <- rbindlist(all_metrics, use.names = TRUE, fill = TRUE)
fwrite(
  combined_metrics,
  file.path(table_dir, "KCNE2_KCNE3_metrics_global_vs_tumor_combined.tsv"),
  sep = "\t"
)

message("Done. KCNE2 and KCNE3 feature plots and metrics exported.")
