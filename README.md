# GSE150290 gastric cancer KCN analysis

This repository contains a single-cell RNA-seq analysis of gastric cancer epithelial states.

The question is simple: are **KCNQ1**, **KCNE2** and **KCNE3** linked to intestinal-like tumor cells, or to diffuse-like tumor cells?

The work was done during a six-month Master 1 bioinformatics internship in the **Regulations of Ion Channel in Cancer** laboratory. The team studies how ion channels contribute to cancer biology, with projects focused on PDAC and gastric cancer. The laboratory is part of Universite Cote d'Azur and affiliated with CNRS and Inserm.

## Dataset

The analysis uses processed single-cell RNA-seq matrices from **GSE150290**.

Source paper:

Kim J. et al. *Single-cell analysis of gastric pre-cancerous and cancer lesions reveals cell lineage diversity and intratumoral heterogeneity*. npj Precision Oncology, 2022. DOI: [10.1038/s41698-022-00251-1](https://doi.org/10.1038/s41698-022-00251-1)

The dataset contains gastric normal, pre-cancerous and cancer lesions. The analysis focuses on the non-immune atlas first, then on the broad tumor compartment.

![GSE150290 overview](GC_dataset1_KCN_intestinal_diffuse_removed-1.png)

## Biological Context

This side project is linked to the team's work on ion channels in gastrointestinal cancers. Ion channels are not only transport proteins. They can also contribute to epithelial organization, differentiation, signaling and tumor progression.

This is relevant for **KCNQ1**, because previous work from Raphael Rapetti-Mauss and colleagues showed that KCNQ1 helps maintain epithelial organization. Loss or inhibition of KCNQ1 can favor beta-catenin redistribution, reduced epithelial integrity and a more proliferative phenotype. This makes it biologically meaningful to ask whether KCNQ1 is more associated with an intestinal-like epithelial tumor state.

Reference: [10.1073/pnas.1702913114](https://doi.org/10.1073/pnas.1702913114)

## Pipeline

![GSE150290 pipeline flowchart](GSE150290_dataset1_pipeline_flowchart.png)

The analysis has two complementary branches.

1. **Markers only**
   Start from the global non-immune t-SNE, select broad tumor cells using marker expression, recluster tumor cells with paper-like Seurat settings, then define intestinal-like and diffuse-like tumor poles.

2. **Markers plus CNV**
   Keep the same logic, then add inferCNV as support for malignant epithelial selection. CNV is not used alone. It is used to refine the tumor compartment already supported by markers.

Final refined tumor map:

- **1,614 refined tumor cells**
- **KCNQ1:** 108 positive cells, 6.69%
- **KCNE2:** 55 positive cells, 3.41%
- **KCNE3:** 341 positive cells, 21.13%

## Main Result

KCNQ1 and KCNE3 are oriented toward the intestinal-like tumor state.

In the refined tumor map, KCNQ1 is positively correlated with the intestinal score and negatively correlated with the diffuse score. KCNE3 shows the same pattern, with a stronger diffuse-negative association. KCNE2 is weaker and less stable.

Main Spearman results:

- **KCNQ1:** rho intestinal = 0.21, rho diffuse = -0.24
- **KCNE3:** rho intestinal = 0.20, rho diffuse = -0.34
- **KCNE2:** weak signal, not clearly interpretable

These results support a transcriptomic association. They do not prove channel activity or causality.

## Folder Structure

```text
gastric-kcnq1-gse150290/
|-- bibliography/              source paper and supplementary material
|-- cache/                     GEO processed matrices
|-- pipeline_markers_only/     first analysis branch
|-- pipeline_markers_plus_cnv/ refined analysis branch with inferCNV support
|-- raw_geo/                   raw GEO files kept for traceability
|-- scripts/                   helper scripts
|-- GSE150290_dataset1_pipeline_flowchart.png
`-- README.md
```

## Key Outputs

- `pipeline_markers_only/outputs/tumor_poles/figures/tumor_all_global_intestinal_diffuse_annotated.pdf`
- `pipeline_markers_only/outputs/kcnq1/figures/KCNQ1_featureplot_global_and_tumor.pdf`
- `pipeline_markers_only/outputs/kcne2_kcne3/figures/KCNE2_featureplot_global_and_tumor.pdf`
- `pipeline_markers_only/outputs/kcne2_kcne3/figures/KCNE3_featureplot_global_and_tumor.pdf`
- `pipeline_markers_plus_cnv/outputs/global_cnv_screen/figures/global_tsne_tumor_markers_and_cnv_like_genes.pdf`
- `pipeline_markers_plus_cnv/outputs/refined_kcnq1_kcne2_kcne3/figures/refined_tumor_KCNQ1_KCNE2_KCNE3_multiplot.pdf`
- `pipeline_markers_plus_cnv/outputs/final/GSE150290_pipeline_signature_sources.xlsx`

## Reproduce

The analysis uses R for Seurat, inferCNV, marker scoring, statistics and figures.

The main scripts are:

- `pipeline_markers_only/scripts/run_global_nonimmune_tumor_extraction.R`
- `pipeline_markers_only/scripts/run_all_global_tumor_cells_exact_paper_pipeline.R`
- `pipeline_markers_only/scripts/annotate_broad_tumor_tsne_intestinal_diffuse.R`
- `pipeline_markers_plus_cnv/scripts/run_infercnv_on_existing_global_tsne.R`
- `pipeline_markers_plus_cnv/scripts/refine_tumor_compartment_markers_plus_cnv.R`
- `pipeline_markers_plus_cnv/scripts/plot_refined_kcnq1_kcne2_kcne3_and_spearman.R`

## Interpretation

This dataset gives a first validation signal for the gastric KCN question.

It is useful because the published paper already contains a tumor-cell structure with intestinal, diffuse/EMT and EmyoT-like programs. The reconstruction is close to the paper logic, but it is not an exact reproduction of the original analysis because the full internal code and detailed annotations were not available.

The most robust interpretation is therefore:

- KCNQ1 and KCNE3 are compatible with an intestinal-like malignant epithelial program.
- KCNE2 is detected, but the biological direction is weaker.
- The result should be interpreted as a single-cell transcriptomic association, not as functional proof.

## References

- Kim J. et al. (2022). *Single-cell analysis of gastric pre-cancerous and cancer lesions reveals cell lineage diversity and intratumoral heterogeneity*. npj Precision Oncology. DOI: [10.1038/s41698-022-00251-1](https://doi.org/10.1038/s41698-022-00251-1)
- Rapetti-Mauss R. et al. (2017). KCNQ1 maintains epithelial organization and contributes to Wnt/beta-catenin regulation. DOI: [10.1073/pnas.1702913114](https://doi.org/10.1073/pnas.1702913114)
