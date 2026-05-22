suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260521)

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_path) || !nzchar(script_path)) {
  script_path <- normalizePath("scripts/01_build_gse165816_foot_skin_seurat.R", winslash = "/", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
count_dir <- file.path(root, "data", "raw", "GSE165816", "counts_csv")
qc_path <- file.path(root, "data", "metadata", "GSE165816_file_qc.tsv")
out_dir <- file.path(root, "results", "seurat_gse165816_foot_skin")
fig_dir <- file.path(root, "figures", "seurat_gse165816_foot_skin")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", root)

cell_type_markers <- list(
  keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10", "KRT6A", "KRT16", "KRT17"),
  fibroblast_stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "RAMP2", "ESAM"),
  myeloid = c("LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "TYROBP", "CTSS"),
  t_nk = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"),
  b_plasma = c("MS4A1", "CD79A", "CD74", "MZB1", "JCHAIN"),
  pericyte_smc = c("RGS5", "PDGFRB", "MCAM", "ACTA2", "TAGLN", "MYH11"),
  melanocyte = c("MLANA", "PMEL", "TYR", "DCT"),
  schwann = c("SOX10", "MPZ", "PLP1", "S100B")
)

program_markers <- list(
  inflammatory_arrest = c("IL1B", "TNF", "CXCL8", "CXCL2", "CCL2", "CCL3", "CCL4", "S100A8", "S100A9", "NFKBIA", "PTGS2"),
  angiogenesis = c("PECAM1", "VWF", "KDR", "FLT1", "ESAM", "EMCN", "ANGPT2", "PLVAP"),
  ecm_remodeling = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "MMP1", "MMP2", "MMP3", "MMP9", "MMP11", "TIMP1"),
  epithelial_migration = c("KRT6A", "KRT6B", "KRT16", "KRT17", "ITGA3", "ITGB1", "LAMC2", "MMP9", "AREG", "HBEGF"),
  hypoxia_oxidative_stress = c("HIF1A", "VEGFA", "SOD2", "HMOX1", "NQO1", "TXN", "JUN", "FOS"),
  senescence_sasp = c("CDKN1A", "CDKN2A", "SERPINE1", "IGFBP7", "MMP3", "CXCL8", "CCL2")
)

read_counts_csv_sparse <- function(path, sample_code) {
  message("Reading ", basename(path))
  barcodes <- strsplit(readLines(gzfile(path), n = 1), ",", fixed = TRUE)[[1]]
  dt <- data.table::fread(path, skip = 1, header = FALSE, data.table = TRUE, showProgress = FALSE)
  genes <- make.unique(dt[[1]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "integer"
  rownames(mat) <- genes
  colnames(mat) <- paste(sample_code, barcodes, sep = "_")
  Matrix(mat, sparse = TRUE)
}

make_sample_object <- function(row) {
  path <- file.path(count_dir, row$file)
  mat <- read_counts_csv_sparse(path, row$sample_code)
  obj <- CreateSeuratObject(
    counts = mat,
    project = "GSE165816_DFU",
    min.cells = 0,
    min.features = 0
  )
  obj$geo_accession <- row$geo_accession
  obj$sample_code <- row$sample_code
  obj$disease <- row$char_disease
  obj$tissue <- row$char_tissue
  obj$sample_title <- row$title
  obj
}

qc <- data.table::fread(qc_path, data.table = FALSE)
foot_qc <- qc %>%
  filter(char_tissue == "Foot skin") %>%
  arrange(match(char_disease, c("Non-diabetic", "Non-DFU Diabetic", "DFU-healer", "DFU-nonhealer")), sample_code)

message("Foot-skin samples: ", nrow(foot_qc))
message("Expected cells from files: ", sum(foot_qc$cells_from_header))

rds_path <- file.path(out_dir, "GSE165816_foot_skin_seurat_firstpass.rds")

if (file.exists(rds_path)) {
  message("Loading existing object: ", rds_path)
  seu <- readRDS(rds_path)
} else {
  objects <- lapply(seq_len(nrow(foot_qc)), function(i) make_sample_object(foot_qc[i, ]))
  message("Merging ", length(objects), " sample objects")
  seu <- Reduce(function(x, y) merge(x, y), objects)
  rm(objects)
  gc()

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$pass_qc <- seu$nFeature_RNA >= 200 & seu$nCount_RNA >= 500 & seu$percent.mt <= 25

  qc_by_sample <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    group_by(geo_accession, sample_code, disease, tissue) %>%
    summarise(
      cells_raw = n(),
      cells_pass_qc = sum(pass_qc),
      median_nFeature_RNA = median(nFeature_RNA),
      median_nCount_RNA = median(nCount_RNA),
      median_percent_mt = median(percent.mt),
      .groups = "drop"
    )
  write.table(qc_by_sample, file.path(out_dir, "GSE165816_foot_skin_qc_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  message("Cells before QC: ", ncol(seu))
  seu <- subset(seu, subset = pass_qc)
  message("Cells after QC: ", ncol(seu))

  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 40, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

  saveRDS(seu, rds_path)
  message("Saved ", rds_path)
}

# Seurat v5 keeps merged samples as multiple assay layers. Join them for
# downstream scoring and plotting.
seu <- JoinLayers(seu, assay = "RNA")

# Manual module scores on log-normalized expression. This avoids AddModuleScore
# failures when a short marker list has too few present genes in a given object.
expr <- GetAssayData(seu, assay = "RNA", layer = "data")
score_gene_set <- function(markers) {
  present <- intersect(markers, rownames(expr))
  if (length(present) == 0) {
    return(rep(0, ncol(expr)))
  }
  Matrix::colMeans(expr[present, , drop = FALSE])
}

marker_presence <- bind_rows(
  lapply(names(cell_type_markers), function(name) {
    tibble::tibble(
      score_group = "cell_type",
      score_name = name,
      requested_genes = paste(cell_type_markers[[name]], collapse = ","),
      present_genes = paste(intersect(cell_type_markers[[name]], rownames(expr)), collapse = ","),
      n_requested = length(cell_type_markers[[name]]),
      n_present = length(intersect(cell_type_markers[[name]], rownames(expr)))
    )
  }),
  lapply(names(program_markers), function(name) {
    tibble::tibble(
      score_group = "program",
      score_name = name,
      requested_genes = paste(program_markers[[name]], collapse = ","),
      present_genes = paste(intersect(program_markers[[name]], rownames(expr)), collapse = ","),
      n_requested = length(program_markers[[name]]),
      n_present = length(intersect(program_markers[[name]], rownames(expr)))
    )
  })
)
write.table(marker_presence, file.path(out_dir, "GSE165816_marker_presence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (name in names(cell_type_markers)) {
  seu[[paste0("ct_", name)]] <- score_gene_set(cell_type_markers[[name]])
}

for (name in names(program_markers)) {
  seu[[paste0("program_", name)]] <- score_gene_set(program_markers[[name]])
}

ct_cols <- paste0("ct_", names(cell_type_markers))
ct_mat <- as.matrix(seu@meta.data[, ct_cols, drop = FALSE])
best <- max.col(ct_mat, ties.method = "first")
seu$broad_cell_type <- names(cell_type_markers)[best]
seu$broad_cell_type_score <- ct_mat[cbind(seq_len(nrow(ct_mat)), best)]
seu$broad_cell_type[seu$broad_cell_type_score < 0] <- "unknown"

saveRDS(seu, rds_path)

meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell")
program_cols <- paste0("program_", names(program_markers))

celltype_counts <- meta %>%
  count(disease, broad_cell_type, name = "cells") %>%
  group_by(disease) %>%
  mutate(fraction_within_disease = cells / sum(cells)) %>%
  ungroup()
write.table(celltype_counts, file.path(out_dir, "GSE165816_foot_skin_celltype_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sample_celltype_program <- meta %>%
  group_by(geo_accession, sample_code, disease, broad_cell_type) %>%
  summarise(
    cells = n(),
    across(all_of(program_cols), mean, .names = "{.col}_mean"),
    .groups = "drop"
  )
write.table(sample_celltype_program, file.path(out_dir, "GSE165816_foot_skin_sample_celltype_program_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

disease_celltype_program <- sample_celltype_program %>%
  group_by(disease, broad_cell_type) %>%
  summarise(
    samples = n_distinct(sample_code),
    cells = sum(cells),
    across(ends_with("_mean"), mean, .names = "{.col}_sample_mean"),
    .groups = "drop"
  )
write.table(disease_celltype_program, file.path(out_dir, "GSE165816_foot_skin_disease_celltype_program_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

compare_healer_nonhealer <- sample_celltype_program %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  tidyr::pivot_longer(cols = ends_with("_mean"), names_to = "program", values_to = "score") %>%
  group_by(broad_cell_type, program) %>%
  summarise(
    n_healer = sum(disease == "DFU-healer"),
    n_nonhealer = sum(disease == "DFU-nonhealer"),
    mean_healer = mean(score[disease == "DFU-healer"], na.rm = TRUE),
    mean_nonhealer = mean(score[disease == "DFU-nonhealer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = mean_healer - mean_nonhealer,
    p_wilcox = ifelse(n_healer >= 2 && n_nonhealer >= 2, wilcox.test(score ~ disease)$p.value, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(p_wilcox)
write.table(compare_healer_nonhealer, file.path(out_dir, "GSE165816_healer_vs_nonhealer_sample_level_wilcox.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

Idents(seu) <- "broad_cell_type"

pdf(file.path(fig_dir, "GSE165816_foot_skin_qc_violin.pdf"), width = 11, height = 4)
print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "disease", pt.size = 0, ncol = 3))
dev.off()

pdf(file.path(fig_dir, "GSE165816_foot_skin_umap_broad_cell_type.pdf"), width = 8, height = 6)
print(DimPlot(seu, reduction = "umap", group.by = "broad_cell_type", raster = TRUE, label = TRUE, repel = TRUE) + NoLegend())
dev.off()

pdf(file.path(fig_dir, "GSE165816_foot_skin_umap_disease.pdf"), width = 8, height = 6)
print(DimPlot(seu, reduction = "umap", group.by = "disease", raster = TRUE))
dev.off()

pdf(file.path(fig_dir, "GSE165816_foot_skin_celltype_fraction_by_disease.pdf"), width = 9, height = 5)
print(
  ggplot(celltype_counts, aes(x = disease, y = fraction_within_disease, fill = broad_cell_type)) +
    geom_col(width = 0.8) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Fraction of QC-passed foot-skin cells", fill = "Broad cell type")
)
dev.off()

focus <- disease_celltype_program %>%
  filter(broad_cell_type %in% c("keratinocyte", "fibroblast_stromal", "endothelial", "myeloid"))

pdf(file.path(fig_dir, "GSE165816_foot_skin_sample_level_program_scores.pdf"), width = 12, height = 8)
for (program in program_cols) {
  p <- sample_celltype_program %>%
    filter(broad_cell_type %in% c("keratinocyte", "fibroblast_stromal", "endothelial", "myeloid")) %>%
    ggplot(aes(x = disease, y = .data[[paste0(program, "_mean")]], color = disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    geom_jitter(width = 0.15, size = 1.8) +
    facet_wrap(~ broad_cell_type, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = program, title = program)
  print(p)
}
dev.off()

sessionInfo_path <- file.path(out_dir, "sessionInfo.txt")
sink(sessionInfo_path)
print(sessionInfo())
sink()

message("Done. Outputs in:")
message("  ", out_dir)
message("  ", fig_dir)
