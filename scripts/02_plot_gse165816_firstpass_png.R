suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_path) || !nzchar(script_path)) {
  script_path <- normalizePath("scripts/02_plot_gse165816_firstpass_png.R", winslash = "/", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "results", "seurat_gse165816_foot_skin")
fig_dir <- file.path(root, "figures", "seurat_gse165816_foot_skin")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(file.path(out_dir, "GSE165816_foot_skin_seurat_firstpass.rds"))
seu <- JoinLayers(seu, assay = "RNA")

png(file.path(fig_dir, "GSE165816_foot_skin_umap_broad_cell_type.png"), width = 1800, height = 1400, res = 180)
print(DimPlot(seu, reduction = "umap", group.by = "broad_cell_type", raster = TRUE, label = TRUE, repel = TRUE) + NoLegend())
dev.off()

png(file.path(fig_dir, "GSE165816_foot_skin_umap_disease.png"), width = 1800, height = 1400, res = 180)
print(DimPlot(seu, reduction = "umap", group.by = "disease", raster = TRUE))
dev.off()

marker_features <- c(
  "KRT14", "KRT5", "KRT10", "KRT16",
  "COL1A1", "COL1A2", "DCN", "LUM",
  "PECAM1", "VWF", "CLDN5",
  "LYZ", "LST1", "S100A8", "S100A9",
  "CD3D", "NKG7",
  "MS4A1", "CD79A", "MZB1",
  "RGS5", "ACTA2", "TAGLN",
  "MLANA", "SOX10"
)
marker_features <- intersect(marker_features, rownames(seu))
png(file.path(fig_dir, "GSE165816_foot_skin_broad_celltype_marker_dotplot.png"), width = 2400, height = 1100, res = 180)
print(
  DotPlot(seu, features = marker_features, group.by = "broad_cell_type") +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 9))
)
dev.off()

sample_scores <- read.delim(file.path(out_dir, "GSE165816_foot_skin_sample_celltype_program_scores.tsv"), check.names = FALSE)
focus_celltypes <- c("keratinocyte", "fibroblast_stromal", "endothelial", "myeloid")
focus_programs <- c(
  "program_inflammatory_arrest_mean",
  "program_ecm_remodeling_mean",
  "program_epithelial_migration_mean",
  "program_hypoxia_oxidative_stress_mean"
)
plot_df <- sample_scores %>%
  filter(broad_cell_type %in% focus_celltypes) %>%
  select(sample_code, disease, broad_cell_type, all_of(focus_programs)) %>%
  tidyr::pivot_longer(cols = all_of(focus_programs), names_to = "program", values_to = "score")

png(file.path(fig_dir, "GSE165816_foot_skin_key_program_scores_by_sample.png"), width = 2400, height = 1600, res = 180)
print(
  ggplot(plot_df, aes(x = disease, y = score, color = disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.15) +
    geom_jitter(width = 0.15, size = 1.6) +
    facet_grid(program ~ broad_cell_type, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = "Mean log-normalized module score per sample/cell type")
)
dev.off()

message("PNG plots written to ", fig_dir)
