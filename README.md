# Stabilizing G-quadruplexes with pyridostatin augments response to immune checkpoint blockade by modulating the tumor immune microenvironment

> **Status:** manuscript under revision (Nature Communications) — if you use this analysis, please cite this repository (see [`CITATION.cff`](CITATION.cff)) until the paper is published.

Kosiol, Kienzl, Gerritzen, Corvino, *et al.* Corresponding author: Annkristin Heine
(Annkristin.Heine@ukbonn.de). Co-senior authors: Katrin Paeschke & Annkristin Heine.

## Summary

Immune checkpoint blockade (ICB) fails in many patients whose tumors are immunologically
"cold." This study shows that stabilizing G-quadruplexes (G4s) — four-stranded DNA secondary
structures — with the G4 ligand **pyridostatin (PDS)** converts cold tumors into inflamed, ICB-
responsive ones. In a checkpoint-resistant B16OVA melanoma model, intratumoral PDS slows tumor
growth, prolongs survival, and sensitizes tumors to anti-PD-1 therapy in a CD8⁺ T cell-dependent
manner. Mechanistically, PDS increases tumor-cell antigen presentation and immunogenic cell
death, drives intratumoral myeloid cells (neutrophils, monocytes, macrophages, DCs) toward a
pro-inflammatory, cytokine/chemokine-producing state, and boosts the recruitment, proliferation,
and cytotoxicity of tumor-infiltrating T and NK cells. Effects on myeloid and tumor cells were
confirmed in human cells in vitro.

## Scope of this repository

This repository contains **only the single-cell RNA-seq analysis** (Dillon Corvino's
contribution: Institute of Experimental Oncology, University Hospital Bonn) — CD45⁺ and CD45⁻
cells sorted from B16OVA tumors 4 days post-treatment. It does **not** contain the bulk 3'
mRNA-seq analyses (tumor cells at day 9, Fig 2a-b; in vitro-activated CD8⁺ T cells at day 3, Fig
4h) or any of the in vivo, flow cytometry, or wet-lab experiments — those live with their
respective authors.

Manuscript panels produced from this repository:

| Panel(s) | Content |
|---|---|
| Fig 5c–d | UMAP of CD45⁺ tumor-infiltrating cells, Vehicle vs. PDS, day 4 |
| Fig 5e | Pseudobulk heatmap of cytokine/chemokine expression across CD45⁺ cell types |
| Fig 5f–g | GO/KEGG enrichment of Monocyte and Neutrophil PDS-vs-Vehicle DEGs |
| Sup Fig 2a–h | ISG and cGAS/STING module-score heatmaps/violins, myeloid vs. CD45⁻ |
| Sup Fig 4b | CD45⁺ cell-subtype composition (day 4) |
| Sup Fig 4c | CD45⁻ cytokine/chemokine pseudobulk heatmap |
| Sup Fig 4d–e | GO/KEGG enrichment of Macrophage and DC PDS-vs-Vehicle DEGs |

See [`figure_manifest.csv`](figure_manifest.csv) for the exact output file, generating script,
and input object behind every one of these panels.

## Repository structure

```
scripts/
├── 01_Nils_dataset_preprocessing.Rmd     # raw CellRanger -> QC/integrated/clustered object
├── 02_CD45Pos_data_annotation.Rmd        # SingleR annotation -> saves/CD45pos_seurat.rds
├── 03_CD45Neg_data_annotation.Rmd        # SingleR + SCEVAN annotation -> saves/CD45neg_seurat.rds
├── manuscript_figures.R                  # regenerates every panel in the table above
└── archive/                              # exploratory analyses & superseded drafts (not
                                           # used by any manuscript panel; kept for reference)
figure_manifest.csv                       # panel -> output file -> script -> input object
```

`saves/` (processed Seurat objects) and `results/` (figures/tables) are not version-controlled —
see [Data availability](#data-availability) below.

## Reproducing the figures

With `saves/CD45pos_seurat.rds` and `saves/CD45neg_seurat.rds` in place:

```bash
Rscript scripts/manuscript_figures.R
```

Outputs are written to `results/manuscript_plots/` as `FigXX_description.{pdf,csv}`.

## Methods (single-cell analysis, as reported in the manuscript)

Ambient RNA was estimated and removed with **SoupX** (v1.6.2). Cells were retained with
200–5,500 detected features, 1,000–45,000 UMIs, <10% mitochondrial reads, and
log10(genes)/UMI > 0.8 (CD45⁺) / > 0.75 (CD45⁻). Data were normalized with **SCTransform**
(`vst.flavor = "v2"`, regressing out mitochondrial content) and integrated across samples by
canonical correlation analysis. Clustering resolution (0.4) was chosen by iterative testing with
**Clustree** (v0.5.1); clusters were annotated by canonical marker expression, and two
basophils identified at day 4 were excluded from day-4-specific analyses.

Differential expression used the Wilcoxon rank-sum test via **Seurat**'s `FindMarkers()`
(v5.3.0). Functional enrichment used **enrichR** (v3.4) and **clusterProfiler** (v4.16) against
the KEGG 2019 Mouse and GO Biological Process 2021 databases, with pathways of interest manually
curated and plotted with **ggplot2**. Pseudobulk expression was computed with
`AggregateExpression()` and visualized with `DoHeatmap()` (**Seurat** v5.5.0), in places using
the Batlow palette via **scico** (v1.5.0); helper/visualization functions from **scCustomize**
(v3.0.1). Per-cell ISG and cGAS/STING module scores used `AddModuleScore()` (Seurat v5.5.0) and
were visualized with `VlnPlot2()` (**SeuratExtend** v1.2.4); Vehicle-vs-PDS comparisons within
each cell type used a two-sided Wilcoxon rank-sum test (`compare_means()`, **ggpubr** v0.6.3)
with Holm correction for multiple comparisons.

## Data availability

CellRanger output and processed Seurat objects are archived on **Zenodo**:
**https://doi.org/10.5281/zenodo.21391578** (restricted access until publication). Raw
sequencing reads will be deposited to ENA. Download the processed data from Zenodo into `data/`
and `saves/` to run the analysis.

> Note: the manuscript's current "Availability of data and materials" section states that
> sequencing data can be provided upon request from the corresponding author, rather than
> naming Zenodo/ENA explicitly. Flagging this so the eventual submitted statement and this
> README can be reconciled before publication.

## Citation

This is unpublished work. Until the paper is published, please cite this repository — see
[`CITATION.cff`](CITATION.cff).

---
Analysis by **Dillon Corvino** · [GitHub](https://github.com/Eomesodermin) · [dilloncorvino.com](https://dilloncorvino.com)
