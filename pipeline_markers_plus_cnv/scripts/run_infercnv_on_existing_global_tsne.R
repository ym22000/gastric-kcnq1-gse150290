suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(infercnv)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Add an inferCNV validation layer on top of the existing gastric pipeline.
# The current t-SNE objects are kept unchanged. We only rebuild the same 13,022
# selected cells with the full raw gene space so that immune markers and CNV
# inference can be used to validate the tumor compartment.

project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
marker_pipeline_dir <- file.path(project_dir, "pipeline_markers_only")
cnv_pipeline_dir <- file.path(project_dir, "pipeline_markers_plus_cnv")

cache_dir <- file.path(project_dir, "cache")
raw_dir <- file.path(project_dir, "raw_geo", "raw_matrices")

global_rds <- file.path(marker_pipeline_dir, "outputs", "global_nonimmune_atlas", "objects", "global_nonimmune_seurat.rds")
tumor_rds <- file.path(
  marker_pipeline_dir,
  "outputs",
  "tumor_reclustering",
  "objects",
  "tumor_all_global_exact_paper_seurat.rds"
)
pole_labels_tsv <- file.path(
  marker_pipeline_dir,
  "outputs",
  "tumor_poles",
  "tables",
  "intestinal_diffuse_cell_labels.tsv"
)
tumor_cells_tsv <- file.path(marker_pipeline_dir, "outputs", "global_nonimmune_atlas", "tables", "all_global_tumor_cells.tsv")

output_dir <- file.path(cnv_pipeline_dir, "outputs", "infercnv_validation")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")
infercnv_dir <- file.path(output_dir, "infercnv_run")
tmp_dir <- file.path(output_dir, "tmp_extract")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(infercnv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

full_gene_rds <- file.path(object_dir, "selected_cells_full_gene_space_seurat.rds")
infercnv_rds <- file.path(object_dir, "infercnv_result_object.rds")

read_processed_cell_ids <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
}

read_raw_selected_cells <- function(raw_tar_path, sample_title, selected_ids) {
  sample_tmp_dir <- file.path(tmp_dir, sample_title)
  dir.create(sample_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  untar(raw_tar_path, exdir = sample_tmp_dir)

  matrix_path <- list.files(sample_tmp_dir, pattern = "matrix\\.mtx$", recursive = TRUE, full.names = TRUE)
  genes_path <- list.files(sample_tmp_dir, pattern = "genes\\.tsv$", recursive = TRUE, full.names = TRUE)
  barcodes_path <- list.files(sample_tmp_dir, pattern = "barcodes\\.tsv$", recursive = TRUE, full.names = TRUE)

  if (length(matrix_path) == 0 || length(genes_path) == 0 || length(barcodes_path) == 0) {
    stop("Could not locate matrix/genes/barcodes files for ", sample_title)
  }

  matrix_path <- matrix_path[1]
  genes_path <- genes_path[1]
  barcodes_path <- barcodes_path[1]

  genes_dt <- fread(genes_path, header = FALSE)
  barcodes <- fread(barcodes_path, header = FALSE)$V1
  raw_mat <- readMM(matrix_path)

  raw_ids <- paste0(sample_title, "-", barcodes)
  keep_idx <- match(selected_ids, raw_ids)
  keep_idx <- keep_idx[!is.na(keep_idx)]
  if (length(keep_idx) != length(selected_ids)) {
    stop("Could not match all selected cell IDs in raw matrix for ", sample_title)
  }

  mat <- raw_mat[, keep_idx, drop = FALSE]
  rownames(mat) <- genes_dt$V2
  colnames(mat) <- selected_ids
  mat
}

collapse_duplicate_genes <- function(mat) {
  gene_names <- rownames(mat)
  if (!anyDuplicated(gene_names)) {
    return(mat)
  }

  message("Keeping one row per duplicated gene symbol for speed...")
  keep <- !duplicated(gene_names)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- gene_names[keep]
  mat
}

make_gene_order_table <- function(gene_symbols) {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_ranges <- genes(txdb)
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = names(gene_ranges),
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )
  gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]
  gene_df <- data.frame(
    ENTREZID = names(gene_ranges),
    chr = as.character(seqnames(gene_ranges)),
    start = start(gene_ranges),
    stop = end(gene_ranges),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, gene_map, by = "ENTREZID")
  gene_df <- gene_df[gene_df$chr %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
  gene_df <- gene_df[gene_df$SYMBOL %in% gene_symbols, ]
  gene_df <- gene_df[order(gene_df$SYMBOL, gene_df$chr, gene_df$start), ]

  # Keep one genomic interval per symbol for inferCNV compatibility.
  gene_df <- gene_df[!duplicated(gene_df$SYMBOL), c("SYMBOL", "chr", "start", "stop")]
  colnames(gene_df) <- c("gene", "chr", "start", "stop")
  gene_df
}

score_cells <- function(obj, markers) {
  markers <- intersect(markers, rownames(obj[["RNA"]]))
  if (length(markers) == 0) {
    return(rep(0, ncol(obj)))
  }
  Matrix::colMeans(GetAssayData(obj, assay = "RNA", layer = "data")[markers, , drop = FALSE])
}

if (!file.exists(full_gene_rds)) {
  message("Rebuilding the selected 13,022 cells with the full raw gene space...")
  current_global <- readRDS(global_rds)

  processed_files <- list.files(cache_dir, pattern = "^GSM.*\\.txt\\.gz$", full.names = TRUE)
  if (length(processed_files) == 0) {
    stop("No processed GSM files found in ", cache_dir)
  }

  mats <- list()
  for (proc_path in processed_files) {
    sample_name <- sub("\\.txt\\.gz$", "", basename(proc_path))
    sample_title <- sub("^[^_]+_", "", sample_name)
    selected_ids <- read_processed_cell_ids(proc_path)
    raw_tar <- file.path(raw_dir, paste0(sample_name, ".raw_gene_bc_matrices.tar.gz"))
    if (!file.exists(raw_tar)) {
      stop("Missing raw matrix tarball: ", raw_tar)
    }
    message("Reading full raw counts for ", sample_title, "...")
    mats[[sample_title]] <- read_raw_selected_cells(raw_tar, sample_title, selected_ids)
  }

  full_counts <- do.call(cbind, mats)
  full_counts <- collapse_duplicate_genes(full_counts)
  full_counts <- full_counts[rowSums(full_counts) > 0, , drop = FALSE]
  full_counts <- full_counts[, colnames(current_global), drop = FALSE]

  full_gene_obj <- CreateSeuratObject(
    counts = full_counts,
    meta.data = current_global@meta.data[colnames(full_counts), , drop = FALSE],
    project = "GSE150290_selected_full_gene_space"
  )
  full_gene_obj <- NormalizeData(full_gene_obj, verbose = FALSE)

  tsne_embeddings <- Embeddings(current_global, "tsne")[colnames(full_gene_obj), , drop = FALSE]
  full_gene_obj[["tsne"]] <- CreateDimReducObject(
    embeddings = tsne_embeddings,
    key = "tSNE_",
    assay = "RNA"
  )

  saveRDS(full_gene_obj, full_gene_rds)
} else {
  full_gene_obj <- readRDS(full_gene_rds)
}

current_global <- readRDS(global_rds)
tumor_obj <- readRDS(tumor_rds)
tumor_cells <- fread(tumor_cells_tsv)$cell_id
pole_labels <- fread(pole_labels_tsv)

full_gene_obj$in_current_tumor_compartment <- colnames(full_gene_obj) %in% tumor_cells
full_gene_obj$tumor_screen_group <- ifelse(
  full_gene_obj$in_current_tumor_compartment,
  "Extracted tumor cells",
  "Other selected cells"
)

monocyte_markers <- c("PTPRC", "LST1", "FCN1", "AIF1", "CTSS", "TYMP", "HLA-DRA", "LYZ", "S100A8", "S100A9", "IL1B")
epithelial_markers <- c("EPCAM", "KRT19", "KRT8", "KRT18", "KRT17", "MUC1", "TFF3", "TFF2", "MUC6")
stromal_markers <- c("COL1A1", "COL3A1", "TAGLN", "DCN", "PDGFRB", "PDGFRA")

full_gene_obj$monocyte_score <- score_cells(full_gene_obj, monocyte_markers)
full_gene_obj$epithelial_score <- score_cells(full_gene_obj, epithelial_markers)
full_gene_obj$stromal_score <- score_cells(full_gene_obj, stromal_markers)

raw_counts <- GetAssayData(full_gene_obj, assay = "RNA", layer = "counts")
positive_marker_count <- rep(0L, ncol(full_gene_obj))
for (gene in intersect(c("PTPRC", "LST1", "FCN1", "AIF1", "CTSS", "TYMP", "HLA-DRA"), rownames(raw_counts))) {
  positive_marker_count <- positive_marker_count + as.integer(raw_counts[gene, ] > 0)
}
full_gene_obj$monocyte_positive_markers <- positive_marker_count

mono_q <- quantile(full_gene_obj$monocyte_score, 0.90, na.rm = TRUE)
epi_q <- quantile(full_gene_obj$epithelial_score, 0.75, na.rm = TRUE)
stromal_q <- quantile(full_gene_obj$stromal_score, 0.75, na.rm = TRUE)

reference_cells <- rownames(full_gene_obj@meta.data)[
  full_gene_obj$monocyte_score >= mono_q &
    full_gene_obj$epithelial_score <= epi_q &
    full_gene_obj$stromal_score <= stromal_q &
    full_gene_obj$monocyte_positive_markers >= 2 &
    !full_gene_obj$in_current_tumor_compartment
]

# Favor adjacent-normal monocyte-like cells when enough are available.
reference_adjacent <- reference_cells[full_gene_obj$tissue_group[reference_cells] == "Adjacent normal"]
if (length(reference_adjacent) >= 50) {
  reference_cells <- reference_adjacent
}

if (length(reference_cells) < 30) {
  stop("Too few monocyte-like reference cells were identified for inferCNV.")
}

full_gene_obj$infercnv_group <- ifelse(
  colnames(full_gene_obj) %in% reference_cells,
  "Monocyte_reference",
  "Observation"
)

ref_summary <- data.frame(
  metric = c(
    "n_selected_cells",
    "n_reference_cells",
    "reference_fraction",
    "reference_adjacent_fraction",
    "monocyte_score_q90_threshold",
    "epithelial_score_q75_threshold",
    "stromal_score_q75_threshold"
  ),
  value = c(
    ncol(full_gene_obj),
    length(reference_cells),
    length(reference_cells) / ncol(full_gene_obj),
    mean(full_gene_obj$tissue_group[reference_cells] == "Adjacent normal"),
    mono_q,
    epi_q,
    stromal_q
  )
)
fwrite(ref_summary, file.path(table_dir, "infercnv_reference_summary.tsv"), sep = "\t")

ref_cells_dt <- data.table(
  cell_id = reference_cells,
  sample_title = full_gene_obj$sample_title[reference_cells],
  tissue_group = full_gene_obj$tissue_group[reference_cells],
  monocyte_score = full_gene_obj$monocyte_score[reference_cells],
  epithelial_score = full_gene_obj$epithelial_score[reference_cells],
  stromal_score = full_gene_obj$stromal_score[reference_cells],
  monocyte_positive_markers = full_gene_obj$monocyte_positive_markers[reference_cells]
)
fwrite(ref_cells_dt, file.path(table_dir, "infercnv_reference_cells.tsv"), sep = "\t")

gene_order_df <- make_gene_order_table(rownames(full_gene_obj[["RNA"]]))
gene_order_df <- gene_order_df[gene_order_df$gene %in% rownames(full_gene_obj[["RNA"]]), ]
gene_order_path <- file.path(table_dir, "hg19_gene_order_for_infercnv.tsv")
fwrite(gene_order_df, gene_order_path, sep = "\t", col.names = FALSE)

annotation_dt <- data.table(
  cell_id = colnames(full_gene_obj),
  group = full_gene_obj$infercnv_group
)
annotation_path <- file.path(table_dir, "infercnv_cell_annotations.tsv")
fwrite(annotation_dt, annotation_path, sep = "\t", col.names = FALSE)

counts_for_infercnv <- raw_counts[gene_order_df$gene, , drop = FALSE]
detected_genes <- rowSums(counts_for_infercnv > 0) >= 5
counts_for_infercnv <- counts_for_infercnv[detected_genes, , drop = FALSE]
gene_order_df <- gene_order_df[gene_order_df$gene %in% rownames(counts_for_infercnv), ]
fwrite(gene_order_df, gene_order_path, sep = "\t", col.names = FALSE)

if (!file.exists(infercnv_rds)) {
  message("Running inferCNV on the selected gastric cells...")
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_for_infercnv,
    annotations_file = annotation_path,
    delim = "\t",
    gene_order_file = gene_order_path,
    ref_group_names = "Monocyte_reference"
  )

  infercnv_res <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = infercnv_dir,
    cluster_by_groups = FALSE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = 8,
    no_plot = TRUE
  )
  saveRDS(infercnv_res, infercnv_rds)
} else {
  infercnv_res <- readRDS(infercnv_rds)
}

expr_cnv <- infercnv_res@expr.data
expr_cnv <- expr_cnv[apply(expr_cnv, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

scale_to_signed_unit <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  -1 + 2 * (x - rng[1]) / diff(rng)
}

expr_scaled <- t(apply(expr_cnv, 1, scale_to_signed_unit))
cnv_score <- sqrt(colMeans(expr_scaled^2, na.rm = TRUE))

cnv_df <- data.table(
  cell_id = names(cnv_score),
  cnv_score = as.numeric(cnv_score)
)
ref_mean <- mean(cnv_df$cnv_score[cnv_df$cell_id %in% reference_cells], na.rm = TRUE)
ref_sd <- sd(cnv_df$cnv_score[cnv_df$cell_id %in% reference_cells], na.rm = TRUE)
cnv_threshold <- ref_mean + 2 * ref_sd
cnv_df[, cnv_high := cnv_score > cnv_threshold]
cnv_df[, infercnv_group := full_gene_obj$infercnv_group[match(cell_id, colnames(full_gene_obj))]]
cnv_df[, sample_title := full_gene_obj$sample_title[match(cell_id, colnames(full_gene_obj))]]
cnv_df[, tissue_group := full_gene_obj$tissue_group[match(cell_id, colnames(full_gene_obj))]]
cnv_df[, in_current_tumor_compartment := cell_id %in% tumor_cells]
fwrite(cnv_df, file.path(table_dir, "infercnv_cell_scores.tsv"), sep = "\t")

tumor_enrichment <- data.table(
  group = c("Extracted tumor cells", "Other selected cells"),
  n_cells = c(
    sum(cnv_df$in_current_tumor_compartment),
    sum(!cnv_df$in_current_tumor_compartment)
  ),
  n_cnv_high = c(
    sum(cnv_df$cnv_high & cnv_df$in_current_tumor_compartment),
    sum(cnv_df$cnv_high & !cnv_df$in_current_tumor_compartment)
  )
)
tumor_enrichment[, pct_cnv_high := 100 * n_cnv_high / n_cells]
fwrite(tumor_enrichment, file.path(table_dir, "infercnv_tumor_enrichment.tsv"), sep = "\t")

full_gene_obj$cnv_score <- cnv_df$cnv_score[match(colnames(full_gene_obj), cnv_df$cell_id)]
full_gene_obj$cnv_high_label <- ifelse(
  cnv_df$cnv_high[match(colnames(full_gene_obj), cnv_df$cell_id)],
  "CNV-high",
  "CNV-low"
)

tumor_obj$cnv_score <- cnv_df$cnv_score[match(colnames(tumor_obj), cnv_df$cell_id)]
tumor_obj$Pole <- pole_labels$PoleLabel[match(colnames(tumor_obj), pole_labels$cell_id)]

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

global_ref_plot <- DimPlot(
  full_gene_obj,
  reduction = "tsne",
  group.by = "infercnv_group",
  cols = c("Observation" = "grey85", "Monocyte_reference" = "#1b9e77"),
  pt.size = 0.45
) +
  labs(
    title = "Global gastric t-SNE",
    subtitle = "Monocyte-like reference cells for inferCNV",
    color = NULL
  ) +
  base_theme

global_tumor_plot <- DimPlot(
  full_gene_obj,
  reduction = "tsne",
  group.by = "tumor_screen_group",
  cols = c("Other selected cells" = "grey85", "Extracted tumor cells" = "#2b6cb0"),
  pt.size = 0.45
) +
  labs(
    title = "Global gastric t-SNE",
    subtitle = "Current extracted tumor compartment",
    color = NULL
  ) +
  base_theme

global_cnv_plot <- FeaturePlot(
  full_gene_obj,
  features = "cnv_score",
  reduction = "tsne",
  cols = c("grey94", "#7b2cbf"),
  pt.size = 0.45
) +
  labs(
    title = "Global gastric t-SNE",
    subtitle = "inferCNV score"
  ) +
  base_theme

global_cnv_high_plot <- DimPlot(
  full_gene_obj,
  reduction = "tsne",
  group.by = "cnv_high_label",
  cols = c("CNV-low" = "grey85", "CNV-high" = "#d73027"),
  pt.size = 0.45
) +
  labs(
    title = "Global gastric t-SNE",
    subtitle = paste0("CNV-high threshold = reference mean + 2 SD (", sprintf("%.3f", cnv_threshold), ")"),
    color = NULL
  ) +
  base_theme

tumor_pole_plot <- DimPlot(
  tumor_obj,
  reduction = "tsne",
  group.by = "Pole",
  cols = c("Diffuse cells" = "#6a4c93", "Intestinal cells" = "#d62828"),
  pt.size = 0.8
) +
  labs(
    title = "Tumor-only gastric t-SNE",
    subtitle = "Final intestinal versus diffuse poles",
    color = NULL
  ) +
  base_theme

tumor_cnv_plot <- FeaturePlot(
  tumor_obj,
  features = "cnv_score",
  reduction = "tsne",
  cols = c("grey94", "#7b2cbf"),
  pt.size = 0.8
) +
  labs(
    title = "Tumor-only gastric t-SNE",
    subtitle = "inferCNV score within the tumor map"
  ) +
  base_theme

pdf(
  file.path(figure_dir, "infercnv_validation_global_and_tumor_tsne.pdf"),
  width = 14,
  height = 7,
  onefile = TRUE
)
print(global_ref_plot + global_tumor_plot)
print(global_cnv_plot + global_cnv_high_plot)
print(tumor_pole_plot + tumor_cnv_plot)
dev.off()

if (dir.exists(tmp_dir)) {
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}

message("Done. inferCNV validation exported.")
