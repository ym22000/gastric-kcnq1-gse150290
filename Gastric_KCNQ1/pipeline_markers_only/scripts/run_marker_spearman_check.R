suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

project_dir <- normalizePath(
  "C:/Users/pcyou/Desktop/stage_LLM/M1_Bioinformatics_Ion_Channel/Gastric_KCNQ1",
  winslash = "/",
  mustWork = TRUE
)
pipeline_dir <- file.path(project_dir, "pipeline_markers_only")

obj <- readRDS(file.path(
  pipeline_dir,
  "outputs",
  "tumor_reclustering",
  "objects",
  "tumor_all_global_exact_paper_seurat.rds"
))

lab <- fread(file.path(
  pipeline_dir,
  "outputs",
  "tumor_poles",
  "tables",
  "intestinal_diffuse_cell_labels.tsv"
))

intestinal_score <- lab[["IntestinalPoleScore"]][match(colnames(obj), lab[["cell_id"]])]
diffuse_score <- lab[["DiffusePoleScore"]][match(colnames(obj), lab[["cell_id"]])]

genes <- c("REG4", "COL1A1")

for (gene in genes) {
  x <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")[gene, ])
  spearman_int <- cor.test(x, intestinal_score, method = "spearman", exact = FALSE)
  spearman_diff <- cor.test(x, diffuse_score, method = "spearman", exact = FALSE)

  cat("\nGENE\t", gene, "\n", sep = "")
  cat(
    "vs_intestinal\trho=", unname(spearman_int$estimate),
    "\tp=", spearman_int$p.value, "\n",
    sep = ""
  )
  cat(
    "vs_diffuse\trho=", unname(spearman_diff$estimate),
    "\tp=", spearman_diff$p.value, "\n",
    sep = ""
  )
}
