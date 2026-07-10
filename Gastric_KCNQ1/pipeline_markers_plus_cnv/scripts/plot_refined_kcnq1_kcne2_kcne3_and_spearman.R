suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
cnv_pipeline_dir <- file.path(project_dir, "pipeline_markers_plus_cnv")

input_rds <- file.path(
  cnv_pipeline_dir,
  "outputs",
  "tumor_refined_markers_cnv",
  "objects",
  "tumor_refined_markers_cnv_seurat.rds"
)

output_dir <- file.path(cnv_pipeline_dir, "outputs", "refined_kcnq1_kcne2_kcne3")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
genes <- c("KCNQ1", "KCNE2", "KCNE3")
missing_genes <- setdiff(genes, rownames(obj))
if (length(missing_genes) > 0) {
  stop("Missing genes in refined object: ", paste(missing_genes, collapse = ", "))
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

pole_plot <- DimPlot(
  obj,
  reduction = "tsne",
  group.by = "PoleLabel",
  cols = c("Intestinal cells" = "#d7301f", "Diffuse cells" = "#7b3294"),
  pt.size = 0.8
) +
  labs(
    title = "Refined tumor t-SNE",
    subtitle = "Intestinal versus diffuse poles",
    color = NULL
  ) +
  base_theme

plot_gene <- function(gene) {
  FeaturePlot(
    obj,
    features = gene,
    reduction = "tsne",
    cols = c("grey94", "#03045e"),
    pt.size = 0.8
  ) +
    labs(
      title = gene,
      subtitle = NULL
    ) +
    base_theme
}

gene_plots <- lapply(genes, plot_gene)

pdf(
  file.path(figure_dir, "refined_tumor_KCNQ1_KCNE2_KCNE3_multiplot.pdf"),
  width = 14,
  height = 10,
  onefile = TRUE
)
print(wrap_plots(c(list(pole_plot), gene_plots), ncol = 2))
dev.off()

results <- list()
for (gene in genes) {
  raw_counts <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "counts")[gene, ])
  norm_values <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")[gene, ])
  positive_raw <- raw_counts[raw_counts > 0]

  spearman_int <- cor.test(norm_values, obj$IntestinalPoleScore, method = "spearman", exact = FALSE)
  spearman_diff <- cor.test(norm_values, obj$DiffusePoleScore, method = "spearman", exact = FALSE)

  results[[length(results) + 1]] <- data.frame(
    gene = gene,
    n_cells = length(raw_counts),
    n_positive_cells = sum(raw_counts > 0),
    pct_positive_cells = 100 * mean(raw_counts > 0),
    total_raw_transcripts = sum(raw_counts),
    mean_raw_counts_per_cell = mean(raw_counts),
    mean_raw_counts_positive_cells = if (length(positive_raw) > 0) mean(positive_raw) else 0,
    max_raw_count = max(raw_counts),
    mean_log_normalized_expression = mean(norm_values),
    spearman_vs_intestinal_rho = unname(spearman_int$estimate),
    spearman_vs_intestinal_p = spearman_int$p.value,
    spearman_vs_diffuse_rho = unname(spearman_diff$estimate),
    spearman_vs_diffuse_p = spearman_diff$p.value,
    stringsAsFactors = FALSE
  )
}

results_df <- rbindlist(results)
fwrite(results_df, file.path(table_dir, "refined_KCNQ1_KCNE2_KCNE3_spearman_metrics.tsv"), sep = "\t")

lines_out <- c("Refined tumor map: KCNQ1, KCNE2, KCNE3",
               "")
for (i in seq_len(nrow(results_df))) {
  row <- results_df[i, ]
  lines_out <- c(
    lines_out,
    paste0(row$gene),
    paste0("positive_cells\t", row$n_positive_cells, "/", row$n_cells, " (", sprintf("%.2f", row$pct_positive_cells), "%)"),
    paste0("total_raw_transcripts\t", row$total_raw_transcripts),
    paste0("mean_raw_counts_per_cell\t", sprintf("%.4f", row$mean_raw_counts_per_cell)),
    paste0("mean_raw_counts_positive_cells\t", sprintf("%.4f", row$mean_raw_counts_positive_cells)),
    paste0("mean_log_normalized_expression\t", sprintf("%.4f", row$mean_log_normalized_expression)),
    paste0("spearman_vs_intestinal\trho=", sprintf("%.4f", row$spearman_vs_intestinal_rho), "\tp=", signif(row$spearman_vs_intestinal_p, 4)),
    paste0("spearman_vs_diffuse\trho=", sprintf("%.4f", row$spearman_vs_diffuse_rho), "\tp=", signif(row$spearman_vs_diffuse_p, 4)),
    ""
  )
}
writeLines(lines_out, file.path(table_dir, "refined_KCNQ1_KCNE2_KCNE3_spearman_metrics.txt"))

message("Done. Refined KCNQ1/KCNE2/KCNE3 projections and Spearman results exported.")
