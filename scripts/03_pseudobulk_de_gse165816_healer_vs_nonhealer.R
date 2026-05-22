suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_path) || !nzchar(script_path)) {
  script_path <- normalizePath("scripts/03_pseudobulk_de_gse165816_healer_vs_nonhealer.R", winslash = "/", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
seurat_path <- file.path(root, "results", "seurat_gse165816_foot_skin", "GSE165816_foot_skin_seurat_firstpass.rds")
out_dir <- file.path(root, "results", "pseudobulk_de_gse165816")
fig_dir <- file.path(root, "figures", "pseudobulk_de_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

focus_celltypes <- c("keratinocyte", "fibroblast_stromal", "endothelial", "myeloid")

message("Loading ", seurat_path)
seu <- readRDS(seurat_path)
seu <- JoinLayers(seu, assay = "RNA")
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
meta <- seu@meta.data %>% tibble::rownames_to_column("cell")

make_pseudobulk <- function(cell_type) {
  sub_meta <- meta %>%
    dplyr::filter(
      broad_cell_type == cell_type,
      disease %in% c("DFU-healer", "DFU-nonhealer")
    )
  sample_counts <- sub_meta %>%
    dplyr::count(sample_code, geo_accession, disease, name = "cells") %>%
    dplyr::filter(cells >= 20)
  sub_meta <- sub_meta %>%
    dplyr::filter(sample_code %in% sample_counts$sample_code)
  samples <- sample_counts$sample_code
  pb <- sapply(samples, function(sample) {
    cell_ids <- sub_meta$cell[sub_meta$sample_code == sample]
    Matrix::rowSums(counts[, cell_ids, drop = FALSE])
  })
  colnames(pb) <- samples
  storage.mode(pb) <- "integer"
  list(counts = pb, sample_meta = sample_counts)
}

run_deseq <- function(cell_type, pb) {
  count_mat <- pb$counts
  sample_meta <- as.data.frame(pb$sample_meta)
  rownames(sample_meta) <- sample_meta$sample_code
  count_mat <- count_mat[, rownames(sample_meta), drop = FALSE]
  keep <- rowSums(count_mat >= 10) >= 3 & rowSums(count_mat) >= 50
  count_mat <- count_mat[keep, , drop = FALSE]
  sample_meta$disease <- factor(sample_meta$disease, levels = c("DFU-nonhealer", "DFU-healer"))
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = sample_meta,
    design = ~ disease
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("disease", "DFU-healer", "DFU-nonhealer"))
  out <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj, pvalue) %>%
    mutate(
      cell_type = cell_type,
      direction = case_when(
        is.na(padj) ~ "NA",
        padj < 0.1 & log2FoldChange > 0 ~ "higher_in_healer",
        padj < 0.1 & log2FoldChange < 0 ~ "higher_in_nonhealer",
        TRUE ~ "not_significant"
      )
    )
  attr(out, "dds") <- dds
  out
}

all_res <- list()
summary_rows <- list()

for (cell_type in focus_celltypes) {
  message("Pseudobulk DE for ", cell_type)
  pb <- make_pseudobulk(cell_type)
  write.table(
    pb$sample_meta,
    file.path(out_dir, paste0("GSE165816_", cell_type, "_pseudobulk_sample_cell_counts.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  if (n_distinct(pb$sample_meta$disease) < 2 || ncol(pb$counts) < 4) {
    warning("Skipping ", cell_type, ": insufficient samples")
    next
  }
  res <- run_deseq(cell_type, pb)
  out_path <- file.path(out_dir, paste0("GSE165816_", cell_type, "_DESeq2_healer_vs_nonhealer.tsv"))
  write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  all_res[[cell_type]] <- res
  summary_rows[[cell_type]] <- tibble::tibble(
    cell_type = cell_type,
    samples_healer = sum(pb$sample_meta$disease == "DFU-healer"),
    samples_nonhealer = sum(pb$sample_meta$disease == "DFU-nonhealer"),
    cells_healer = sum(pb$sample_meta$cells[pb$sample_meta$disease == "DFU-healer"]),
    cells_nonhealer = sum(pb$sample_meta$cells[pb$sample_meta$disease == "DFU-nonhealer"]),
    genes_tested = nrow(res),
    sig_padj_0_1_higher_healer = sum(res$padj < 0.1 & res$log2FoldChange > 0, na.rm = TRUE),
    sig_padj_0_1_higher_nonhealer = sum(res$padj < 0.1 & res$log2FoldChange < 0, na.rm = TRUE)
  )

  plot_df <- res %>%
    mutate(
      neglog10padj = -log10(pmax(padj, 1e-300)),
      hit = padj < 0.1 & abs(log2FoldChange) > 0.5
    )
  png(file.path(fig_dir, paste0("GSE165816_", cell_type, "_volcano_healer_vs_nonhealer.png")), width = 1200, height = 1000, res = 160)
  print(
    ggplot(plot_df, aes(x = log2FoldChange, y = neglog10padj, color = hit)) +
      geom_point(size = 0.7, alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#b2182b")) +
      theme_classic() +
      labs(
        title = paste0(cell_type, ": DFU-healer vs DFU-nonhealer"),
        x = "log2FC (healer / non-healer)",
        y = "-log10 adjusted P"
      ) +
      theme(legend.position = "none")
  )
  dev.off()
}

combined <- bind_rows(all_res)
write.table(combined, file.path(out_dir, "GSE165816_all_focus_celltypes_DESeq2_healer_vs_nonhealer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- bind_rows(summary_rows)
write.table(summary, file.path(out_dir, "GSE165816_pseudobulk_DE_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_hits <- combined %>%
  dplyr::filter(!is.na(padj)) %>%
  group_by(cell_type) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 50) %>%
  ungroup()
write.table(top_hits, file.path(out_dir, "GSE165816_top50_DE_genes_by_focus_celltype.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Done. Outputs in ", out_dir)
print(summary)
