suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Step 1 of the final pipeline:
# build the global non-immune atlas from the processed GSM matrices,
# assign broad lineages with simple marker programs,
# and highlight the global tumor-annotated cells used downstream.
project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
pipeline_dir <- file.path(project_dir, "pipeline_markers_only")
cache_dir <- file.path(project_dir, "cache")
output_dir <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

sample_files <- list.files(cache_dir, pattern = "^GSM.*\\.txt\\.gz$", full.names = TRUE)
if (length(sample_files) == 0) {
  stop("No GSM matrices found in ", cache_dir)
}

read_sample_matrix <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  dt <- data.table::fread(path, sep = "\t", header = FALSE, skip = 1, data.table = FALSE)
  genes <- dt[[1]]
  expr <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- genes
  colnames(expr) <- header
  Matrix(expr, sparse = TRUE)
}

build_cell_metadata <- function(cell_ids, sample_title, geo_accession) {
  # The local GSM files correspond to the non-immune cells retained in the paper.
  # Using the same convention as the published subset: -A are matched adjacent sites,
  # -B are gastric cancer lesions.
  tissue_group <- if (grepl("-A$", sample_title)) "Adjacent normal" else "Cancer"
  data.frame(
    cell_id = cell_ids,
    sample_title = sample_title,
    geo_accession = geo_accession,
    tissue_group = tissue_group,
    stringsAsFactors = FALSE,
    row.names = cell_ids
  )
}

run_paper_like_seurat <- function(seu) {
  seu <- NormalizeData(seu, verbose = FALSE)
  # Closest Seurat v5 equivalent of the paper's Seurat mean/dispersion filters:
  # x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 0.5
  seu <- FindVariableFeatures(
    seu,
    selection.method = "mean.var.plot",
    mean.cutoff = c(0.0125, 6),
    dispersion.cutoff = c(0.5, Inf),
    verbose = FALSE
  )
  # The supplementary describes scaling with UMI regression for selected populations.
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

score_marker_program <- function(avg_mat, genes) {
  genes <- intersect(genes, rownames(avg_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(avg_mat)))
  }
  colMeans(avg_mat[genes, , drop = FALSE])
}

message("Reading processed GSE150290 non-immune matrices...")
matrix_list <- list()
meta_list <- list()
reference_genes <- NULL

for (path in sample_files) {
  sample_name <- sub("\\.txt\\.gz$", "", basename(path))
  geo_accession <- sub("_.*$", "", sample_name)
  sample_title <- sub("^[^_]+_", "", sample_name)
  mat <- read_sample_matrix(path)

  if (is.null(reference_genes)) {
    reference_genes <- rownames(mat)
  } else if (!identical(reference_genes, rownames(mat))) {
    stop("Gene order differs in ", basename(path))
  }

  matrix_list[[sample_title]] <- mat
  meta_list[[sample_title]] <- build_cell_metadata(colnames(mat), sample_title, geo_accession)
}

combined_counts <- do.call(cbind, matrix_list)
combined_meta <- bind_rows(meta_list)
combined_meta <- combined_meta[colnames(combined_counts), , drop = FALSE]

message("Building global non-immune Seurat object...")
global_obj <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = combined_meta,
  project = "GSE150290_nonimmune_global"
)
global_obj <- run_paper_like_seurat(global_obj)

write.table(
  data.frame(
    metric = c("n_cells", "n_genes", "adjacent_cells", "cancer_cells"),
    value = c(
      ncol(global_obj),
      nrow(global_obj),
      sum(global_obj$tissue_group == "Adjacent normal"),
      sum(global_obj$tissue_group == "Cancer")
    )
  ),
  file.path(table_dir, "global_nonimmune_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

marker_programs <- list(
  tumor = c("EPCAM", "CDH17", "COL3A1", "PDGFRB"),
  IM = c("TFF3", "CDX1", "CDX2"),
  PMC = c("GKN1", "GKN2", "MUC5AC"),
  GMC = c("MUC6", "TFF2"),
  fibroblast = c("MMP2", "PDGFRA", "MYL9", "FN1", "CAV1"),
  EC = c("PLVAP", "KDR", "PTPRB"),
  enteroendocrine = c("CHGA", "GAST", "PROX1")
)

features_for_summary <- intersect(unique(unlist(marker_programs)), rownames(global_obj))
avg_expr <- AverageExpression(
  global_obj,
  features = features_for_summary,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)$RNA
avg_expr <- as.data.frame(avg_expr)
cluster_cols <- colnames(avg_expr)

cluster_summary <- global_obj@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    adjacent_fraction = mean(tissue_group == "Adjacent normal"),
    cancer_fraction = mean(tissue_group == "Cancer"),
    .groups = "drop"
  )

for (program_name in names(marker_programs)) {
  cluster_summary[[paste0(program_name, "_score")]] <- score_marker_program(avg_expr, marker_programs[[program_name]])[
    match(paste0("g", cluster_summary$cluster), cluster_cols)
  ]
}

# Each global cluster receives the lineage with the highest average marker score.
score_cols <- paste0(names(marker_programs), "_score")
cluster_summary$predicted_lineage <- names(marker_programs)[max.col(cluster_summary[, score_cols], ties.method = "first")]

marker_columns_to_keep <- c("EPCAM", "CDH17", "COL3A1", "PDGFRB", "TFF3", "CDX1", "CDX2", "GKN1", "GKN2", "MUC5AC", "MUC6", "TFF2", "MMP2", "PDGFRA", "PLVAP", "KDR", "PTPRB", "CHGA", "GAST", "PROX1")
for (gene in intersect(marker_columns_to_keep, rownames(avg_expr))) {
  cluster_summary[[gene]] <- as.numeric(avg_expr[gene, paste0("g", cluster_summary$cluster)])
}

# Conservative tumor compartment without CNV: published tumor-marker program must dominate,
# and the cluster must be enriched in cancer lesions.
cluster_summary$is_tumor_compartment <- with(
  cluster_summary,
  predicted_lineage == "tumor" & cancer_fraction >= 0.50
)

cluster_summary <- cluster_summary %>%
  arrange(desc(is_tumor_compartment), desc(cancer_fraction), desc(tumor_score))

write.table(
  cluster_summary,
  file.path(table_dir, "global_cluster_marker_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

broad_tumor_clusters <- as.character(cluster_summary$cluster[cluster_summary$predicted_lineage == "tumor"])
if (length(broad_tumor_clusters) == 0) {
  stop("No globally annotated tumor clusters were identified.")
}

broad_tumor_cells <- rownames(global_obj@meta.data)[as.character(global_obj$seurat_clusters) %in% broad_tumor_clusters]

write.table(
  data.frame(
    cell_id = broad_tumor_cells,
    cluster = as.character(global_obj$seurat_clusters[broad_tumor_cells]),
    sample_title = global_obj$sample_title[broad_tumor_cells],
    tissue_group = global_obj$tissue_group[broad_tumor_cells],
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "all_global_tumor_cells.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(global_tumor_clusters = broad_tumor_clusters),
  file.path(table_dir, "selected_global_tumor_clusters.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    metric = c("n_global_tumor_clusters", "global_tumor_clusters", "n_global_tumor_cells"),
    value = c(length(broad_tumor_clusters), paste(broad_tumor_clusters, collapse = ","), length(broad_tumor_cells))
  ),
  file.path(table_dir, "all_global_tumor_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

global_cluster_plot <- DimPlot(
  global_obj,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.45
) +
  labs(
    title = "GSE150290 non-immune global t-SNE",
    subtitle = "Paper-like Seurat parameters | clusters",
    color = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

lineage_map <- setNames(cluster_summary$predicted_lineage, cluster_summary$cluster)
global_obj$predicted_lineage <- factor(
  unname(lineage_map[as.character(global_obj$seurat_clusters)]),
  levels = names(marker_programs)
)

global_lineage_plot <- DimPlot(
  global_obj,
  reduction = "tsne",
  group.by = "predicted_lineage",
  label = FALSE,
  pt.size = 0.45
) +
  labs(
    title = "GSE150290 non-immune global t-SNE",
    subtitle = "Marker-based lineage labels from the paper",
    color = "Lineage"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

global_tumor_highlight <- DimPlot(
  global_obj,
  reduction = "tsne",
  cells.highlight = broad_tumor_cells,
  cols = "grey85",
  cols.highlight = "#2b6cb0",
  pt.size = 0.45
) +
  labs(
    title = "Global tumor-annotated cells",
    subtitle = paste0("Broad tumor set extracted from global t-SNE: ", length(broad_tumor_cells), " cells"),
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figure_dir, "global_nonimmune_tsne_and_tumor_extraction.pdf"),
  global_cluster_plot + global_lineage_plot + plot_layout(ncol = 1),
  width = 9.4,
  height = 10.4,
  device = cairo_pdf
)

ggsave(
  file.path(figure_dir, "global_nonimmune_tumor_cells_highlight_blue.pdf"),
  global_tumor_highlight,
  width = 8.2,
  height = 6.8,
  device = cairo_pdf
)

saveRDS(global_obj, file.path(object_dir, "global_nonimmune_seurat.rds"))

message("Done. Global clustering and tumor compartment extraction completed.")
