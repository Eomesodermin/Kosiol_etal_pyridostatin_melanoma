# scripts/manuscript_figures.R
#
# Self-contained script that regenerates exactly the manuscript-figure panels
# for the pyridostatin B16OVA scRNA-seq analysis: Fig 5c-g, Sup Fig 2a-h,
# Sup Fig 4b-e. Code consolidated (verbatim logic, no styling/parameter
# changes) from scripts/02, 03, 04, 05. See figure_manifest.csv at the repo
# root for full per-panel provenance (new script + archived original).
#
# INPUT -- this script LOADS the already-annotated Seurat objects; it does
# NOT regenerate them. Those objects are produced by two upstream stages,
# both intentionally left untouched by this cleanup:
#   Stage 1 (raw -> integrated/clustered/QC'd):
#     scripts/01_Nils_dataset_preprocessing.Rmd
#     -> saves/Preprocessed_CD45Pos.rds, saves/Preprocessed_CD45Neg.rds
#   Stage 2 (SingleR/SCEVAN annotation, cluster_anno assignment):
#     scripts/02_CD45Pos_data_annotation.Rmd, scripts/03_CD45Neg_data_annotation.Rmd
#     -> saves/CD45pos_seurat.rds, saves/CD45neg_seurat.rds   <-- loaded below
#
# OUTPUT -- results/manuscript_plots/FigXX_description.{pdf,csv}
#
# Run from anywhere inside the repo: Rscript scripts/manuscript_figures.R
# (working directory is auto-detected from the script's own path -- no
# hardcoded absolute paths, portable for GitHub/Zenodo).

## ---- portable working-directory bootstrap ----
.get_script_dir <- function() {
  cmd.args <- commandArgs(trailingOnly = FALSE)
  file.flag <- "--file="
  match.idx <- grep(file.flag, cmd.args)
  if (length(match.idx) > 0) {
    return(dirname(normalizePath(sub(file.flag, "", cmd.args[match.idx]))))
  }
  getwd()  # fallback for interactive/sourced use: assumes cwd is repo root
}
repo.root <- dirname(.get_script_dir())
setwd(repo.root)

set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(scCustomize)
  library(usefulfunctions)
  library(scico)
  library(enrichR)
  library(scales)
  library(SeuratExtend)
  library(ggpubr)
})

output.dir <- "results/manuscript_plots/"
if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

# Headless equivalent of the original chunks' interactive
# `print(p); dev.copy(pdf, path); dev.off()` pattern (needs a live RStudio
# device). Same plot object, same content.
save_pdf <- function(plot_expr, path, width = 10, height = 8) {
  pdf(path, width = width, height = height)
  print(plot_expr)
  dev.off()
}

cat(sprintf("[%s] manuscript_figures.R starting. repo.root = %s\n", Sys.time(), repo.root))


########################################################################
## Fig 5c, 5d, 5e, 5f, 5g  +  Sup Fig 4b, 4d, 4e
## Source: scripts/04_CD45Pos_Manuscript_plots.Rmd
## Input:  saves/CD45pos_seurat.rds
########################################################################

cat(sprintf("[%s] Loading CD45pos.seurat...\n", Sys.time()))
CD45pos.seurat <- LoadSeuratRds("saves/CD45pos_seurat.rds")
DefaultAssay(CD45pos.seurat) <- "SCT"
CD45pos.seurat@meta.data$Treatment <- factor(CD45pos.seurat@meta.data$Treatment, levels = c("Vehicle", "PDS"))

# Day4-only subset, basophils removed (only 2 basophils at Day4 -- too few
# for a meaningful plot/comparison; matches 04's "set_up_object" chunk).
Idents(CD45pos.seurat) <- CD45pos.seurat@meta.data$Time
day4.seurat <- subset(CD45pos.seurat, idents = "Day4")
Idents(day4.seurat) <- day4.seurat@meta.data$cell_type
day4.seurat <- subset(day4.seurat, idents = "Basophils", invert = TRUE)


## ---- Fig 5c: UMAP, Vehicle, Day4 (downsampled 7000/treatment for fair visual comparison) ----
Idents(day4.seurat) <- day4.seurat@meta.data$Treatment
small.seurat <- subset(day4.seurat, downsample = 7000)

Idents(small.seurat) <- small.seurat@meta.data$Treatment
temp.seurat <- subset(small.seurat, idents = "Vehicle")
p.veh <- scCustomize::DimPlot_scCustom(temp.seurat,
  reduction = "umap", label = TRUE, figure_plot = TRUE,
  group.by = "cell_type", pt.size = 1)
save_pdf(p.veh, paste0(output.dir, "Fig5c_UMAP_vehicle.pdf"))

## ---- Fig 5d: UMAP, PDS, Day4 ----
Idents(small.seurat) <- small.seurat@meta.data$Treatment
temp.seurat <- subset(small.seurat, idents = "PDS")
p.pds <- scCustomize::DimPlot_scCustom(temp.seurat,
  reduction = "umap", label = TRUE, figure_plot = TRUE,
  group.by = "cell_type", pt.size = 1)
save_pdf(p.pds, paste0(output.dir, "Fig5d_UMAP_PDS.pdf"))
rm(small.seurat, temp.seurat)


## ---- Sup Fig 4b: cluster-frequency CSV (panel itself is made in Prism from this CSV -- no PDF needed) ----
dist.table <- table(day4.seurat@meta.data$cell_type, day4.seurat@meta.data$Treatment)
treatment.total <- colSums(dist.table)
freq.table <- sweep(dist.table, 2, treatment.total, FUN = "/")
freq.table <- freq.table * 100
freq.table <- as.data.frame(freq.table)
freq.table$Var2 <- factor(freq.table$Var2, levels = c("Vehicle", "PDS"))
write.csv(freq.table, paste0(output.dir, "SupFig4b_cluster_freq_day4.csv"))


## ---- Fig 5e: pseudobulk heatmap, custom cytokine/chemokine panel, Day4, cell_type x Treatment ----
DefaultAssay(day4.seurat) <- "SCT"
pseudo.cd45pos <- AggregateExpression(day4.seurat, assays = "SCT", return.seurat = TRUE,
  group.by = c("cell_type", "Treatment"))
custom.genes <- c("Il1b", "Il6", "Il12a", "Il12b", "Ifng", "Tnf", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Cxcl1", "Cxcl9", "Cxcl10")
# Cxcl13 not present in CD45Pos dataset (per original 04 comment)

p.h1 <- DoHeatmap(pseudo.cd45pos, features = custom.genes, draw.lines = FALSE, angle = 90, raster = FALSE) +
  ggplot2::scale_fill_gradient2(low = "blue", high = "red2")
save_pdf(p.h1, paste0(output.dir, "Fig5e_heatmap_custom_genes.pdf"))
rm(pseudo.cd45pos)


## ---- Fig 5f, 5g, Sup Fig 4d, 4e: Day4 myeloid PDS-vs-Vehicle DEGs -> GO/KEGG dot plots ----
## Restricted to the 4 cell types that actually feed a manuscript panel
## (original 04 script computed DEGs for every Day4 cell type present,
## including B cells/NK cells/Stromal cells/T cells, none of which are used
## in any of the 10 target panels here).

day4.seurat <- PrepSCTFindMarkers(day4.seurat)
DefaultAssay(day4.seurat) <- "SCT"
day4.seurat@meta.data$Treatment_cluster <- paste0(day4.seurat@meta.data$Treatment, "_", day4.seurat@meta.data$cell_type)
Idents(day4.seurat) <- "Treatment_cluster"

celltype.vec <- c("Neutrophils", "Monocytes", "Macrophages", "DC")
DEG.list <- list()
for (i in seq_along(celltype.vec)) {
  cat(sprintf("[%s] FindMarkers Day4 %s (PDS vs Vehicle)...\n", Sys.time(), celltype.vec[i]))
  DEG.list[[celltype.vec[i]]] <- FindMarkers(day4.seurat, assay = "SCT",
    ident.1 = paste0("PDS_", celltype.vec[i]),
    ident.2 = paste0("Vehicle_", celltype.vec[i]),
    only.pos = FALSE, recorrect_umi = FALSE, verbose = TRUE)
}

enrichR::setEnrichrSite(site = "Enrichr")
dbs_use <- c("GO_Biological_Process_2021", "KEGG_2019_Mouse")

filter_terms <- function(data, go_terms, group_name, prefix = "") {
  data %>%
    dplyr::select(c(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score)) %>%
    filter(grepl(paste(go_terms, collapse = "|"), Term, ignore.case = TRUE)) %>%
    mutate(Term = paste0(prefix, Term), group = group_name)
}

process_data <- function(input.1, input.2, input.3, input.4, pval_threshold = 0.05) {
  combined.term.table <- bind_rows(input.1, input.2, input.3, input.4)
  combined.term.table <- combined.term.table %>%
    mutate(Overlap_Ratio = as.numeric(sub("/.*", "", Overlap)) / as.numeric(sub(".*/", "", Overlap))) %>%
    mutate(color_value = ifelse(Adjusted.P.value > pval_threshold, NA, Combined.Score)) %>%
    mutate(color_value_log = log1p(color_value))
  combined.term.table$group <- factor(combined.term.table$group, levels = c("PDS", "Vehicle"))
  combined.term.table <- combined.term.table %>%
    arrange(group, desc(Combined.Score)) %>%
    mutate(Term_order = factor(Term, levels = unique(Term)))
  combined.term.table
}

# NOTE: original 04 also had a plot_all() (unfiltered "select_GOs.pdf")
# variant; that output isn't one of the 10 manuscript panels, so it's
# omitted here (Step 3: "regenerate every panel ... and only those panels").
plot_pds_sig <- function(term.table, title.var = "", pval_threshold = 0.05, save.path, width.var = 12, height.var = 8) {
  significant_pds_terms <- term.table %>%
    filter(group == "PDS" & Adjusted.P.value < pval_threshold) %>%
    pull(Term)
  filtered_combined_table <- term.table %>% filter(Term %in% significant_pds_terms)
  p1 <- ggplot(filtered_combined_table, aes(x = group, y = Term_order, size = Overlap_Ratio)) +
    geom_point(aes(color = color_value_log), na.rm = TRUE) +
    scale_color_gradient(low = "blue", high = "red", na.value = "gray",
      name = "Combined Score", labels = scales::label_log()) +
    scale_size_continuous(range = c(0, 10)) +
    scale_y_discrete(limits = rev) +
    theme_minimal() +
    labs(title = title.var, x = "Treatment group (Day 4)", y = NULL)
  ggsave(save.path, p1, width = width.var, height = height.var)
}

run_celltype_go_dotplot <- function(celltype, GO.terms, KEGG.terms, out.filename) {
  cat(sprintf("[%s] GO/KEGG dotplot: %s\n", Sys.time(), celltype))
  up.genes <- DEG.list[[celltype]] %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(avg_log2FC > 0) %>% rownames()
  dn.genes <- DEG.list[[celltype]] %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(avg_log2FC < 0) %>% rownames()

  up.enriched.terms <- enrichR::enrichr(up.genes, dbs_use)
  dn.enriched.terms <- enrichR::enrichr(dn.genes, dbs_use)

  GO.terms.rx <- paste0("^", GO.terms, " ")
  KEGG.terms.rx <- paste0("^", KEGG.terms, "$")

  PDS.go <- filter_terms(up.enriched.terms$GO_Biological_Process_2021, GO.terms.rx, "PDS")
  Vehcile.go <- filter_terms(dn.enriched.terms$GO_Biological_Process_2021, GO.terms.rx, "Vehicle")
  PDS.kegg <- filter_terms(up.enriched.terms$KEGG_2019_Mouse, KEGG.terms.rx, "PDS", prefix = "KEGG: ")
  Vehcile.kegg <- filter_terms(dn.enriched.terms$KEGG_2019_Mouse, KEGG.terms.rx, "Vehicle", prefix = "KEGG: ")

  term.table.out <- process_data(PDS.go, PDS.kegg, Vehcile.go, Vehcile.kegg)

  plot_pds_sig(term.table.out, title.var = paste0(celltype, " PDS vs. Vehicle"),
    save.path = paste0(output.dir, out.filename))
}

## ---- Fig 5g: Neutrophils GO/KEGG dot plot ----
run_celltype_go_dotplot("Neutrophils",
  GO.terms = c("cytokine-mediated signaling pathway", "neutrophil degranulation", "neutrophil mediated immunity",
    "cellular response to type I interferon", "phagosome acidification", "phagosome maturation",
    "translation", "gene expression", "ribosome biogenesis"),
  KEGG.terms = c("Lysosome", "Phagosome", "NOD-like receptor signaling pathway", "Apoptosis",
    "NF-kappa B signaling pathway", "TNF signaling pathway", "Antigen processing and presentation",
    "RIG-I-like receptor signaling pathway", "Ribosome", "Oxidative phosphorylation"),
  out.filename = "Fig5g_Neutrophils_GO_KEGG_dotplot.pdf")

## ---- Fig 5f: Monocytes GO/KEGG dot plot ----
run_celltype_go_dotplot("Monocytes",
  GO.terms = c("cytokine-mediated signaling pathway", "inflammatory response", "positive regulation of inflammatory response",
    "cellular response to type I interferon", "positive regulation of cytokine production", "translation",
    "gene expression", "ribosome biogenesis"),
  KEGG.terms = c("Lysosome", "Phagosome", "NOD-like receptor signaling pathway", "Apoptosis",
    "Toll-like receptor signaling pathway", "TNF signaling pathway", "Ribosome", "Oxidative phosphorylation"),
  out.filename = "Fig5f_Monocytes_GO_KEGG_dotplot.pdf")

## ---- Sup Fig 4d: Macrophages GO/KEGG dot plot ----
run_celltype_go_dotplot("Macrophages",
  GO.terms = c("cytokine-mediated signaling pathway", "regulation of inflammatory response", "inflammatory response",
    "toll-like receptor signaling pathway", "pattern recognition receptor signaling pathway",
    "cellular response to Interleukin-1", "translation", "gene expression", "ribosome biogenesis"),
  KEGG.terms = c("Lysosome", "Phagosome", "NOD-like receptor signaling pathway", "Apoptosis",
    "Toll-like receptor signaling pathway", "TNF signaling pathway", "Antigen processing and presentation",
    "Ribosome", "Oxidative phosphorylation"),
  out.filename = "SupFig4d_Macrophages_GO_KEGG_dotplot.pdf")

## ---- Sup Fig 4e: DC GO/KEGG dot plot ----
run_celltype_go_dotplot("DC",
  GO.terms = c("cytokine-mediated signaling pathway", "inflammatory response", "cellular response to Interleukin-1",
    "positive regulation of inflammatory response", "positive regulation of cytokine production",
    "ribosome biogenesis", "translation", "gene expression"),
  KEGG.terms = c("Phagosome", "NOD-like receptor signaling pathway", "Apoptosis", "Toll-like receptor signaling pathway",
    "TNF signaling pathway", "NF-kappa B signaling pathway", "RIG-I-like receptor signaling pathway",
    "Ribosome", "Oxidative phosphorylation"),
  out.filename = "SupFig4e_DC_GO_KEGG_dotplot.pdf")

# Keep DEG.list around for Step-4 verification; drop the rest.
CD45pos.DEG.list.day4 <- DEG.list
rm(CD45pos.seurat, day4.seurat, DEG.list)
gc()


########################################################################
## Sup Fig 4c
## Source: scripts/03_CD45Neg_data_analysis.Rmd, chunk chemokine_cytokine_expression
## Input:  saves/CD45neg_seurat.rds
########################################################################

cat(sprintf("[%s] Loading CD45neg.seurat...\n", Sys.time()))
CD45neg.seurat <- LoadSeuratRds("saves/CD45neg_seurat.rds")
DefaultAssay(CD45neg.seurat) <- "SCT"

pseudo.cd45neg <- AggregateExpression(CD45neg.seurat, assays = "SCT", return.seurat = TRUE,
  group.by = c("cluster_anno"))

custom.genes.neg <- c("Il1b", "Il6", "Il12a", "Ifng", "Tnf", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Cxcl1", "Cxcl9", "Cxcl10", "Cxcl13")
# Il12b not present in CD45Neg dataset (per original 03 comment)

p.sup4c <- DoHeatmap(pseudo.cd45neg, features = custom.genes.neg, draw.lines = FALSE, angle = 90, raster = FALSE) +
  scico::scale_fill_scico(palette = "batlow", direction = 1, na.value = "white")
save_pdf(p.sup4c, paste0(output.dir, "SupFig4c_chemokine_heatmap_CD45neg.pdf"))

rm(CD45neg.seurat, pseudo.cd45neg)
gc()


########################################################################
## Sup Fig 2a-h
## Source: scripts/05_Supplementary_Figure_2_ISG_cGAS.Rmd (copied verbatim --
## this was already a complete, self-contained analysis; not re-derived)
## Input:  saves/CD45pos_seurat.rds, saves/CD45neg_seurat.rds
########################################################################

cat(sprintf("[%s] Sup Fig 2: loading objects...\n", Sys.time()))
cd45pos.seurat <- LoadSeuratRds("saves/CD45pos_seurat.rds")
cd45neg.seurat <- LoadSeuratRds("saves/CD45neg_seurat.rds")

## ---- Gene signatures ----
isg.genes <- c(
  "Isg15", "Ifit1", "Ifit2", "Ifit3", "Irf7", "Stat1", "Stat2", "Oas1a",
  "Mx2", "Rsad2", "Usp18", "Ifi44", "Ifi44l", "Oasl1", "Gbp2", "Cxcl10",
  "Ifih1", "Ifi27", "Ifi27l2b", "Ifi30", "Ifi35", "Ifitm1", "Ifitm2", "Ifitm3",
  "Isg20", "Ddx60", "Eif2ak2", "Ly6e", "Oas1g", "Parp12", "Parp14", "Rtp4",
  "Oas1b", "Oas2", "Oas3", "Oasl2", "Ddx58"
)
cgas.genes <- c(
  "Cgas", "Ifi204", "Ddx41", "Tmem173", "Tbk1", "Ikbke", "Irf3", "Irf7",
  "Nfkb1", "Rela", "Trex1", "Enpp1", "Ifnb1", "Il6", "Tnf", "Cxcl10",
  "Ccl2", "Ccl5", "Aim2", "Pycard", "Casp1"
)

## ---- Prepare CD45+ Myeloid Day4 subset ----
myeloid.keep <- c("DC", "Macrophages", "Monocytes", "Neutrophils")
Idents(cd45pos.seurat) <- cd45pos.seurat@meta.data$Time
cd45pos.d4 <- subset(cd45pos.seurat, idents = "Day4")
Idents(cd45pos.d4) <- cd45pos.d4@meta.data$cell_type
cd45pos.d4 <- subset(cd45pos.d4, idents = myeloid.keep)
DefaultAssay(cd45pos.d4) <- "SCT"
cd45pos.d4$Treatment <- factor(cd45pos.d4$Treatment, levels = c("Vehicle", "PDS"))
cd45pos.d4$cell_type <- factor(cd45pos.d4$cell_type, levels = myeloid.keep)
myeloid.col.order <- c("DC_Vehicle", "Macrophages_Vehicle", "Monocytes_Vehicle", "Neutrophils_Vehicle",
                       "DC_PDS", "Macrophages_PDS", "Monocytes_PDS", "Neutrophils_PDS")

## ---- Prepare CD45- Tumor+Fibroblast subset ----
DefaultAssay(cd45neg.seurat) <- "SCT"
cd45neg.seurat$cell_type <- ifelse(
  grepl("^Tumor", cd45neg.seurat$cluster_anno), "Tumor", as.character(cd45neg.seurat$cluster_anno)
)
cd45neg.d4 <- subset(cd45neg.seurat, subset = cell_type %in% c("Tumor", "Fibroblasts"))
cd45neg.d4$Treatment <- factor(cd45neg.d4$Treatment, levels = c("Vehicle", "PDS"))
cd45neg.d4$cell_type <- factor(cd45neg.d4$cell_type, levels = c("Tumor", "Fibroblasts"))
cd45neg.col.order <- c("Tumor_Vehicle", "Fibroblasts_Vehicle", "Tumor_PDS", "Fibroblasts_PDS")

## ---- Helper: ordered pseudobulk heatmap (blue-white-red2, matches Fig5e style) ----
plot_pseudobulk_heatmap <- function(seu, features, group.by, col.order, title) {
  pseudo <- AggregateExpression(seu, assays = "SCT", return.seurat = TRUE, group.by = group.by)
  features.present <- intersect(features, rownames(pseudo))
  col.order <- col.order[col.order %in% colnames(pseudo)]
  Idents(pseudo) <- factor(Idents(pseudo), levels = col.order)
  pseudo <- pseudo[, col.order]
  p <- DoHeatmap(pseudo, features = features.present, draw.lines = FALSE, angle = 90, raster = FALSE) +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red2") +
    ggplot2::labs(title = title)
  list(plot = p, n.features = length(features.present), n.total = length(features))
}

## ---- Helper: violin plot with Holm-adjusted significance stars ----
plot_violin_holm <- function(seu, feature, group.by, split.by, ylab, title,
                              step.increase = 0.12, tip.length = 0.03) {
  p <- VlnPlot2(seu, features = feature, group.by = group.by, split.by = split.by,
    stat.method = "none", pt = FALSE, hide.outlier = TRUE) +
    labs(title = title, y = ylab, x = "Cell Type")

  scores <- p$data
  group_by_arg <- if ("f2" %in% colnames(scores)) c("feature", "f") else "feature"
  formula <- if ("f2" %in% colnames(scores)) value ~ f2 else value ~ f

  stat.test <- ggpubr::compare_means(formula, data = scores, method = "wilcox.test",
    group.by = group_by_arg, p.adjust.method = "holm")

  stat.test$p.adj.signif <- as.character(symnum(stat.test$p.adj,
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
  stat.test <- dplyr::filter(stat.test, p.adj.signif != "ns")

  if (nrow(stat.test) == 0) {
    message("No comparisons significant after Holm adjustment for: ", title)
    return(p)
  }

  stat.test$groups <- apply(stat.test, 1, function(x) c(x[["group1"]], x[["group2"]]), simplify = FALSE)
  stat.test <- SeuratExtend:::vlnplot2_Stat_add_y(stat.test, scores = scores, step.increase = step.increase)
  stat.test$annotation <- stat.test$p.adj.signif

  p + ggpubr::stat_pvalue_manual(stat.test, label = "annotation", tip.length = tip.length)
}

## ---- Panel a: ISG heatmap, Myeloid ----
res <- plot_pseudobulk_heatmap(cd45pos.d4, isg.genes, c("cell_type", "Treatment"), myeloid.col.order, "ISG Signature - Myeloid Cells")
pdf(paste0(output.dir, "SupFig2a_ISG_heatmap_myeloid.pdf"), width = 8, height = 9); print(res$plot); dev.off()

## ---- Panel b: ISG heatmap, CD45- (Tumor + Fibroblasts) ----
res <- plot_pseudobulk_heatmap(cd45neg.d4, isg.genes, c("cell_type", "Treatment"), cd45neg.col.order, "ISG Signature - Tumor & Fibroblasts")
pdf(paste0(output.dir, "SupFig2b_ISG_heatmap_CD45neg.pdf"), width = 6, height = 9); print(res$plot); dev.off()

## ---- Panel c: ISG module score violin, Myeloid ----
isg.present <- intersect(isg.genes, rownames(cd45pos.d4))
cd45pos.d4 <- AddModuleScore(cd45pos.d4, features = list(isg.present), name = "ISG", assay = "SCT")
cd45pos.d4$ISG_score <- cd45pos.d4$ISG1
cd45pos.d4$ISG1 <- NULL
p.isg.vln.myeloid <- plot_violin_holm(cd45pos.d4, "ISG_score", "cell_type", "Treatment", "ISG Score", "ISG Module Score - Myeloid Cells")
ggsave(paste0(output.dir, "SupFig2c_ISG_violin_myeloid.pdf"), p.isg.vln.myeloid, width = 10, height = 6, units = "in")

## ---- Panel d: ISG module score violin, CD45- ----
isg.present.cd45neg <- intersect(isg.genes, rownames(cd45neg.d4))
cd45neg.d4 <- AddModuleScore(cd45neg.d4, features = list(isg.present.cd45neg), name = "ISG", assay = "SCT")
cd45neg.d4$ISG_score <- cd45neg.d4$ISG1
cd45neg.d4$ISG1 <- NULL
p.isg.vln.cd45neg <- plot_violin_holm(cd45neg.d4, "ISG_score", "cell_type", "Treatment", "ISG Score", "ISG Module Score - Tumor & Fibroblasts")
ggsave(paste0(output.dir, "SupFig2d_ISG_violin_CD45neg.pdf"), p.isg.vln.cd45neg, width = 8, height = 6, units = "in")

## ---- Panel e: cGAS/STING heatmap, Myeloid ----
res <- plot_pseudobulk_heatmap(cd45pos.d4, cgas.genes, c("cell_type", "Treatment"), myeloid.col.order, "cGAS/STING Signature - Myeloid Cells")
pdf(paste0(output.dir, "SupFig2e_cGAS_heatmap_myeloid.pdf"), width = 8, height = 7); print(res$plot); dev.off()

## ---- Panel f: cGAS/STING heatmap, CD45- ----
res <- plot_pseudobulk_heatmap(cd45neg.d4, cgas.genes, c("cell_type", "Treatment"), cd45neg.col.order, "cGAS/STING Signature - Tumor & Fibroblasts")
pdf(paste0(output.dir, "SupFig2f_cGAS_heatmap_CD45neg.pdf"), width = 6, height = 7); print(res$plot); dev.off()

## ---- Panel g: cGAS/STING module score violin, Myeloid ----
cgas.present <- intersect(cgas.genes, rownames(cd45pos.d4))
cd45pos.d4 <- AddModuleScore(cd45pos.d4, features = list(cgas.present), name = "cGAS", assay = "SCT")
cd45pos.d4$cGAS_score <- cd45pos.d4$cGAS1
cd45pos.d4$cGAS1 <- NULL
p.cgas.vln.myeloid <- plot_violin_holm(cd45pos.d4, "cGAS_score", "cell_type", "Treatment", "cGAS/STING Score", "cGAS/STING Module Score - Myeloid Cells")
ggsave(paste0(output.dir, "SupFig2g_cGAS_violin_myeloid.pdf"), p.cgas.vln.myeloid, width = 10, height = 6, units = "in")

## ---- Panel h: cGAS/STING module score violin, CD45- ----
cgas.present.cd45neg <- intersect(cgas.genes, rownames(cd45neg.d4))
cd45neg.d4 <- AddModuleScore(cd45neg.d4, features = list(cgas.present.cd45neg), name = "cGAS", assay = "SCT")
cd45neg.d4$cGAS_score <- cd45neg.d4$cGAS1
cd45neg.d4$cGAS1 <- NULL
p.cgas.vln.cd45neg <- plot_violin_holm(cd45neg.d4, "cGAS_score", "cell_type", "Treatment", "cGAS/STING Score", "cGAS/STING Module Score - Tumor & Fibroblasts")
ggsave(paste0(output.dir, "SupFig2h_cGAS_violin_CD45neg.pdf"), p.cgas.vln.cd45neg, width = 8, height = 6, units = "in")

rm(cd45pos.seurat, cd45neg.seurat, cd45pos.d4, cd45neg.d4)
gc()

cat(sprintf("[%s] DONE manuscript_figures.R -- all panels written to %s\n", Sys.time(), output.dir))
