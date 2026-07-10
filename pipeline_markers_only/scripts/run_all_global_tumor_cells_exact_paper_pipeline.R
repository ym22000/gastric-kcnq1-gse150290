suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

# Step 2 of the final pipeline:
# take all globally tumor-annotated cells, rebuild a tumor-only object,
# and run the paper-like Seurat workflow on this broad tumor set.
project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
pipeline_dir <- file.path(project_dir, "pipeline_markers_only")

global_rds <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas", "objects", "global_nonimmune_seurat.rds")
cluster_summary_tsv <- file.path(pipeline_dir, "outputs", "global_nonimmune_atlas", "tables", "global_cluster_marker_summary.tsv")

output_dir <- file.path(pipeline_dir, "outputs", "tumor_reclustering")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

global_obj <- readRDS(global_rds)
cluster_summary <- fread(cluster_summary_tsv)

tumor_clusters <- cluster_summary %>%
  filter(predicted_lineage == "tumor") %>%
  pull(cluster) %>%
  as.character()

if (length(tumor_clusters) == 0) {
  stop("No global clusters are annotated as tumor.")
}

tumor_cells <- rownames(global_obj@meta.data)[as.character(global_obj$seurat_clusters) %in% tumor_clusters]
tumor_counts <- GetAssayData(global_obj, assay = "RNA", layer = "counts")[, tumor_cells, drop = FALSE]
tumor_meta <- global_obj@meta.data[tumor_cells, , drop = FALSE]

obj <- CreateSeuratObject(
  counts = tumor_counts,
  meta.data = tumor_meta,
  project = "GSE150290_all_global_tumor_exact"
)

message("Running exact paper-like tumor-only workflow on all global tumor-annotated cells...")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(
  obj,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.0125, 6),
  dispersion.cutoff = c(0.5, Inf),
  verbose = FALSE
)
obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  vars.to.regress = "nCount_RNA",
  verbose = FALSE
)
obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 30,
  seed.use = 12345,
  verbose = FALSE
)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
obj <- RunTSNE(
  obj,
  dims = 1:5,
  seed.use = 12345,
  check_duplicates = FALSE,
  verbose = FALSE
)

paper_markers <- intersect(
  c("CDH17", "REG4", "MUC13", "IGFBP5", "COL1A1", "S100A4", "TAGLN", "EGR1"),
  rownames(obj)
)

avg_expr <- AverageExpression(
  obj,
  features = paper_markers,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)$RNA
avg_expr <- as.data.frame(avg_expr)
avg_expr$gene <- rownames(avg_expr)

cluster_sizes <- as.data.frame(table(obj$seurat_clusters), stringsAsFactors = FALSE)
colnames(cluster_sizes) <- c("cluster", "n_cells")
cluster_sizes <- cluster_sizes[order(as.integer(cluster_sizes$cluster)), , drop = FALSE]

cluster_marker_summary <- cluster_sizes
for (gene in paper_markers) {
  cluster_marker_summary[[gene]] <- as.numeric(avg_expr[avg_expr$gene == gene, paste0("g", cluster_marker_summary$cluster)])
}

sample_summary <- obj@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  count(cluster, sample_title, tissue_group, name = "n_cells") %>%
  arrange(as.integer(cluster), desc(n_cells))

fwrite(cluster_sizes, file.path(table_dir, "tumor_all_global_cluster_sizes.tsv"), sep = "\t")
fwrite(cluster_marker_summary, file.path(table_dir, "tumor_all_global_cluster_marker_summary.tsv"), sep = "\t")
fwrite(sample_summary, file.path(table_dir, "tumor_all_global_cluster_sample_summary.tsv"), sep = "\t")
write.table(
  data.frame(
    metric = c("n_cells", "n_global_tumor_clusters", "global_tumor_clusters", "n_reclustered_clusters", "variable_features"),
    value = c(ncol(obj), length(tumor_clusters), paste(tumor_clusters, collapse = ","), nlevels(obj$seurat_clusters), length(VariableFeatures(obj)))
  ),
  file.path(table_dir, "tumor_all_global_exact_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(obj, file.path(object_dir, "tumor_all_global_exact_paper_seurat.rds"))

message("Done. Broad tumor-only object and summary tables exported.")
