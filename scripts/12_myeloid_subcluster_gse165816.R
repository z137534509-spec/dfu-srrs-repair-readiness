options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260521)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
seurat_file <- file.path(
  project_dir,
  "results", "seurat_gse165816_foot_skin",
  "GSE165816_foot_skin_seurat_firstpass.rds"
)
out_dir <- file.path(project_dir, "results", "myeloid_subcluster_gse165816")
fig_dir <- file.path(project_dir, "figures", "myeloid_subcluster_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading Seurat object...")
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

message("Subsetting myeloid cells...")
my <- subset(obj, subset = broad_cell_type == "myeloid")
message("Myeloid cells: ", ncol(my))

myeloid_marker_sets <- list(
  classical_monocyte = c("FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "LYZ", "LST1", "CTSS"),
  inflammatory_mono_mac = c("IL1B", "TNF", "CXCL8", "CXCL2", "CCL3", "CCL4", "NFKBIA", "PTGS2"),
  c1q_apoe_macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "APOC1", "MSR1", "MAFB", "MRC1"),
  resolution_metabolic = c("FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CES1", "CXCL5", "TNC", "COL4A1", "COL4A2"),
  antigen_presentation_dc = c("HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "FCER1A", "CLEC10A", "CD1C", "IRF8"),
  langerhans_like = c("CD1A", "CD207", "CD1E", "EPCAM", "FCER1A", "HLA-DQB1"),
  interferon_response = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "IFI6", "IRF7", "OAS1"),
  proliferating = c("MKI67", "TOP2A", "UBE2C", "STMN1", "HMGB2", "TYMS"),
  mast_contamination = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
  neutrophil_like = c("CSF3R", "FCGR3B", "CXCR2", "S100A8", "S100A9", "MMP9", "OLR1")
)

get_data <- function(so) {
  tryCatch(
    GetAssayData(so, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(so, assay = "RNA", slot = "data")
  )
}

score_set <- function(so, genes) {
  mat <- get_data(so)
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) {
    return(rep(0, ncol(so)))
  }
  as.numeric(Matrix::colMeans(mat[present, , drop = FALSE]))
}

message("Running myeloid-only Seurat workflow...")
my <- NormalizeData(my, verbose = FALSE)
my <- FindVariableFeatures(my, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
my <- ScaleData(my, features = VariableFeatures(my), verbose = FALSE)
my <- RunPCA(my, features = VariableFeatures(my), npcs = 30, verbose = FALSE)
my <- FindNeighbors(my, dims = 1:20, verbose = FALSE)
my <- FindClusters(my, resolution = 0.5, verbose = FALSE)
my <- RunUMAP(my, dims = 1:20, verbose = FALSE)
my$myeloid_subcluster <- paste0("M", as.character(Idents(my)))

message("Scoring myeloid states...")
for (nm in names(myeloid_marker_sets)) {
  my[[paste0("myeloid_score_", nm)]] <- score_set(my, myeloid_marker_sets[[nm]])
}

score_cols <- paste0("myeloid_score_", names(myeloid_marker_sets))
cluster_scores <- my@meta.data %>%
  group_by(myeloid_subcluster) %>%
  summarise(
    cells = n(),
    across(all_of(score_cols), mean),
    .groups = "drop"
  )
cluster_scores$assigned_state <- names(myeloid_marker_sets)[
  max.col(as.matrix(cluster_scores[, score_cols]), ties.method = "first")
]
cluster_scores$assigned_state_score <- apply(as.matrix(cluster_scores[, score_cols]), 1, max)
cluster_scores <- cluster_scores %>% arrange(myeloid_subcluster)

state_map <- setNames(cluster_scores$assigned_state, cluster_scores$myeloid_subcluster)
my$myeloid_state <- unname(state_map[my$myeloid_subcluster])

message("Writing summaries...")
saveRDS(my, file.path(out_dir, "GSE165816_myeloid_subcluster_seurat.rds"))

write.table(
  cluster_scores,
  file.path(out_dir, "GSE165816_myeloid_state_marker_scores_by_subcluster.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

counts_by_disease <- my@meta.data %>%
  count(disease, myeloid_subcluster, myeloid_state, name = "cells") %>%
  group_by(disease) %>%
  mutate(fraction_within_disease_myeloid = cells / sum(cells)) %>%
  ungroup()
write.table(
  counts_by_disease,
  file.path(out_dir, "GSE165816_myeloid_subcluster_counts_by_disease.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sample_cluster_counts <- my@meta.data %>%
  count(geo_accession, sample_code, disease, myeloid_subcluster, myeloid_state, name = "cells") %>%
  group_by(geo_accession, sample_code, disease) %>%
  mutate(total_myeloid_cells = sum(cells), fraction_within_sample_myeloid = cells / total_myeloid_cells) %>%
  ungroup()
write.table(
  sample_cluster_counts,
  file.path(out_dir, "GSE165816_myeloid_subcluster_fraction_by_sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sample_state_scores <- my@meta.data %>%
  group_by(geo_accession, sample_code, disease, myeloid_state) %>%
  summarise(
    cells = n(),
    across(all_of(score_cols), mean),
    .groups = "drop"
  )
write.table(
  sample_state_scores,
  file.path(out_dir, "GSE165816_myeloid_state_scores_by_sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sample_all_scores <- my@meta.data %>%
  group_by(geo_accession, sample_code, disease) %>%
  summarise(
    cells = n(),
    across(all_of(score_cols), mean),
    .groups = "drop"
  )
write.table(
  sample_all_scores,
  file.path(out_dir, "GSE165816_myeloid_all_state_scores_by_sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

wilcox_rows <- list()
dfu_samples <- sample_all_scores %>% filter(disease %in% c("DFU-healer", "DFU-nonhealer"))
for (col in score_cols) {
  h <- dfu_samples[[col]][dfu_samples$disease == "DFU-healer"]
  n <- dfu_samples[[col]][dfu_samples$disease == "DFU-nonhealer"]
  wilcox_rows[[col]] <- data.frame(
    feature = col,
    n_healer = length(h),
    n_nonhealer = length(n),
    mean_healer = mean(h),
    mean_nonhealer = mean(n),
    delta_healer_minus_nonhealer = mean(h) - mean(n),
    p_wilcox = suppressWarnings(wilcox.test(h, n, exact = FALSE)$p.value)
  )
}

dfu_cluster <- sample_cluster_counts %>% filter(disease %in% c("DFU-healer", "DFU-nonhealer"))
for (cl in sort(unique(dfu_cluster$myeloid_subcluster))) {
  z <- dfu_cluster %>% filter(myeloid_subcluster == cl)
  h <- z$fraction_within_sample_myeloid[z$disease == "DFU-healer"]
  n <- z$fraction_within_sample_myeloid[z$disease == "DFU-nonhealer"]
  if (length(h) > 0 && length(n) > 0) {
    wilcox_rows[[paste0("fraction_", cl)]] <- data.frame(
      feature = paste0("fraction_", cl),
      n_healer = length(h),
      n_nonhealer = length(n),
      mean_healer = mean(h),
      mean_nonhealer = mean(n),
      delta_healer_minus_nonhealer = mean(h) - mean(n),
      p_wilcox = suppressWarnings(wilcox.test(h, n, exact = FALSE)$p.value)
    )
  }
}
wilcox_summary <- bind_rows(wilcox_rows)
wilcox_summary$padj_bh <- p.adjust(wilcox_summary$p_wilcox, method = "BH")
wilcox_summary <- wilcox_summary %>% arrange(p_wilcox)
write.table(
  wilcox_summary,
  file.path(out_dir, "GSE165816_myeloid_healer_vs_nonhealer_wilcox.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Finding myeloid subcluster markers...")
markers <- FindAllMarkers(
  my,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.25,
  verbose = FALSE
)
markers <- markers %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
write.table(
  markers,
  file.path(out_dir, "GSE165816_myeloid_subcluster_markers.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_head(n = 15) %>%
  ungroup()
write.table(
  top_markers,
  file.path(out_dir, "GSE165816_myeloid_subcluster_top15_markers.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Saving figures...")
png(file.path(fig_dir, "GSE165816_myeloid_umap_subcluster.png"), width = 1800, height = 1500, res = 220)
print(DimPlot(my, reduction = "umap", group.by = "myeloid_subcluster", label = TRUE, repel = TRUE) + NoLegend())
dev.off()

png(file.path(fig_dir, "GSE165816_myeloid_umap_state.png"), width = 1800, height = 1500, res = 220)
print(DimPlot(my, reduction = "umap", group.by = "myeloid_state", label = TRUE, repel = TRUE) + NoLegend())
dev.off()

png(file.path(fig_dir, "GSE165816_myeloid_umap_disease.png"), width = 1800, height = 1500, res = 220)
print(DimPlot(my, reduction = "umap", group.by = "disease"))
dev.off()

marker_features <- unique(c(
  "LYZ", "LST1", "S100A8", "S100A9", "IL1B", "CXCL8",
  "C1QA", "C1QB", "APOE", "CD74", "HLA-DRA", "FCER1A", "CD1C",
  "FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CES1", "CXCL5", "TNC",
  "ISG15", "IFIT1", "MKI67", "TOP2A"
))
marker_features <- intersect(marker_features, rownames(my))
png(file.path(fig_dir, "GSE165816_myeloid_marker_dotplot.png"), width = 2600, height = 1400, res = 220)
print(DotPlot(my, features = marker_features, group.by = "myeloid_subcluster") + RotatedAxis())
dev.off()

feature_genes <- intersect(c("FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CXCL5", "TNC"), rownames(my))
if (length(feature_genes) > 0) {
  png(file.path(fig_dir, "GSE165816_myeloid_resolution_featureplots.png"), width = 2600, height = 1800, res = 220)
  print(FeaturePlot(my, features = feature_genes, reduction = "umap", ncol = 4, order = TRUE))
  dev.off()
}

state_plot <- sample_state_scores %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  mutate(disease = factor(disease, levels = c("DFU-nonhealer", "DFU-healer")))
if (nrow(state_plot) > 0) {
  png(file.path(fig_dir, "GSE165816_myeloid_resolution_score_by_state_sample.png"), width = 2400, height = 1600, res = 220)
  print(
    ggplot(state_plot, aes(x = disease, y = myeloid_score_resolution_metabolic, fill = disease)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
      geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.4, alpha = 0.8) +
      facet_wrap(~ myeloid_state, scales = "free_y") +
      scale_fill_manual(values = c(`DFU-nonhealer` = "#D55E00", `DFU-healer` = "#0072B2")) +
      labs(x = NULL, y = "Resolution-metabolic score", fill = NULL) +
      theme_classic(base_size = 10) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
  )
  dev.off()
}

cluster_frac_plot <- sample_cluster_counts %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  mutate(disease = factor(disease, levels = c("DFU-nonhealer", "DFU-healer")))
if (nrow(cluster_frac_plot) > 0) {
  png(file.path(fig_dir, "GSE165816_myeloid_subcluster_fraction_by_sample.png"), width = 2600, height = 1700, res = 220)
  print(
    ggplot(cluster_frac_plot, aes(x = disease, y = fraction_within_sample_myeloid, fill = disease)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
      geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.4, alpha = 0.8) +
      facet_wrap(~ myeloid_subcluster, scales = "free_y") +
      scale_fill_manual(values = c(`DFU-nonhealer` = "#D55E00", `DFU-healer` = "#0072B2")) +
      labs(x = NULL, y = "Fraction within sample myeloid", fill = NULL) +
      theme_classic(base_size = 10) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
  )
  dev.off()
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done.")
