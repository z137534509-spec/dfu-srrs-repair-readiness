options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
gse_dir <- file.path(project_dir, "data", "raw", "GSE241132")
zip_dir <- file.path(gse_dir, "raw_zip")
tenx_dir <- file.path(gse_dir, "tenx")
meta_file <- file.path(gse_dir, "GSE241132_cell_metadata.txt.gz")

fib_dfu_file <- file.path(project_dir, "results", "fibroblast_trajectory_gse165816", "GSE165816_fibroblast_trajectory_seurat.rds")
out_dir <- file.path(project_dir, "results", "human_wound_roadmap_gse241132")
fig_dir <- file.path(project_dir, "figures", "human_wound_roadmap_gse241132")
dir.create(tenx_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(gzfile(meta_file), check.names = FALSE)
fib_meta <- meta %>%
  filter(newMainCellTypes == "Fibroblast", Doublet == "Singlet") %>%
  mutate(
    timepoint = case_when(
      Condition == "Skin" ~ "D0_skin",
      Condition == "Wound1" ~ "D1_wound",
      Condition == "Wound7" ~ "D7_wound",
      Condition == "Wound30" ~ "D30_wound",
      TRUE ~ Condition
    ),
    timepoint = factor(timepoint, levels = c("D0_skin", "D1_wound", "D7_wound", "D30_wound"))
  )

sample_ids <- sort(unique(fib_meta$orig.ident))

read_sample_fibroblasts <- function(sample_id) {
  sample_dir <- file.path(tenx_dir, sample_id)
  matrix_file <- file.path(sample_dir, paste0(sample_id, "_matrix.mtx.gz"))
  if (!file.exists(matrix_file)) {
    zip_file <- list.files(zip_dir, pattern = paste0(sample_id, "\\.zip$"), full.names = TRUE)
    if (length(zip_file) != 1) {
      stop("Cannot find zip for ", sample_id)
    }
    unzip(zip_file, exdir = tenx_dir)
  }
  mat <- ReadMtx(
    mtx = file.path(sample_dir, paste0(sample_id, "_matrix.mtx.gz")),
    cells = file.path(sample_dir, paste0(sample_id, "_barcodes.tsv.gz")),
    features = file.path(sample_dir, paste0(sample_id, "_features.tsv.gz")),
    feature.column = 2,
    unique.features = TRUE
  )
  colnames(mat) <- paste(sample_id, colnames(mat), sep = "_")
  keep <- intersect(colnames(mat), fib_meta$barcode[fib_meta$orig.ident == sample_id])
  mat[, keep, drop = FALSE]
}

message("Reading GSE241132 fibroblast matrices...")
ref_mats <- lapply(sample_ids, read_sample_fibroblasts)
common_genes <- Reduce(intersect, lapply(ref_mats, rownames))
ref_mats <- lapply(ref_mats, function(x) x[common_genes, , drop = FALSE])
ref_counts <- do.call(cbind, ref_mats)
ref_meta <- fib_meta %>%
  filter(barcode %in% colnames(ref_counts)) %>%
  arrange(match(barcode, colnames(ref_counts)))
stopifnot(identical(ref_meta$barcode, colnames(ref_counts)))
rownames(ref_meta) <- ref_meta$barcode

ref_obj <- CreateSeuratObject(ref_counts, meta.data = ref_meta, project = "GSE241132_fibroblast")
ref_obj <- NormalizeData(ref_obj, verbose = FALSE)
ref_log <- LayerData(ref_obj, assay = "RNA", layer = "data")
ref_counts_layer <- LayerData(ref_obj, assay = "RNA", layer = "counts")

conditions <- levels(ref_meta$timepoint)
avg_log <- sapply(conditions, function(tp) {
  cells <- ref_meta$barcode[ref_meta$timepoint == tp]
  Matrix::rowMeans(ref_log[, cells, drop = FALSE])
})
colnames(avg_log) <- conditions

pct_expr <- sapply(conditions, function(tp) {
  cells <- ref_meta$barcode[ref_meta$timepoint == tp]
  Matrix::rowMeans(ref_counts_layer[, cells, drop = FALSE] > 0)
})
colnames(pct_expr) <- conditions

bad_gene <- grepl("^MT-|^RPL|^RPS|^MALAT1$|^HBB|^HBA", rownames(avg_log), ignore.case = FALSE)

signature_rows <- bind_rows(lapply(conditions, function(tp) {
  other <- setdiff(conditions, tp)
  specificity <- avg_log[, tp] - apply(avg_log[, other, drop = FALSE], 1, max)
  data.frame(
    signature = paste0("GSE241132_", tp, "_fibroblast_specific"),
    gene = rownames(avg_log),
    condition = tp,
    avg_log_expr = avg_log[, tp],
    max_other_avg_log_expr = apply(avg_log[, other, drop = FALSE], 1, max),
    specificity = specificity,
    pct_expr = pct_expr[, tp]
  ) %>%
    filter(!bad_gene, pct_expr >= 0.05, specificity > 0) %>%
    arrange(desc(specificity)) %>%
    slice_head(n = 100)
}))

wound_vs_skin_rows <- bind_rows(lapply(setdiff(conditions, "D0_skin"), function(tp) {
  delta <- avg_log[, tp] - avg_log[, "D0_skin"]
  data.frame(
    signature = paste0("GSE241132_", tp, "_vs_D0_skin_fibroblast_up"),
    gene = rownames(avg_log),
    condition = tp,
    avg_log_expr = avg_log[, tp],
    skin_avg_log_expr = avg_log[, "D0_skin"],
    delta_vs_skin = delta,
    pct_expr = pct_expr[, tp]
  ) %>%
    filter(!bad_gene, pct_expr >= 0.05, delta_vs_skin > 0) %>%
    arrange(desc(delta_vs_skin)) %>%
    slice_head(n = 100)
}))

signature_sets <- split(signature_rows$gene, signature_rows$signature)
signature_sets <- c(signature_sets, split(wound_vs_skin_rows$gene, wound_vs_skin_rows$signature))

fib_dfu <- readRDS(fib_dfu_file)
DefaultAssay(fib_dfu) <- "RNA"
dfu_log <- LayerData(fib_dfu, assay = "RNA", layer = "data")

score_gene_set <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < 3) {
    return(rep(NA_real_, ncol(expr)))
  }
  mat <- as.matrix(expr[genes, , drop = FALSE])
  keep <- apply(mat, 1, sd, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 3) {
    return(rep(NA_real_, ncol(expr)))
  }
  z <- t(scale(t(mat)))
  colMeans(z, na.rm = TRUE)
}

dfu_scores <- as.data.frame(sapply(signature_sets, function(gs) score_gene_set(dfu_log, gs)))
dfu_scores$cell <- rownames(dfu_scores)
dfu_scores <- dfu_scores %>%
  left_join(
    fib_dfu@meta.data %>%
      mutate(cell = rownames(.)) %>%
      select(cell, sample_code, disease, fibro_cluster, repair_axis, fibro_pseudotime),
    by = "cell"
  )

score_cols <- names(signature_sets)
dfu_sample_scores <- dfu_scores %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  mutate(healing_status = ifelse(disease == "DFU-healer", "Healer", "Non-healer")) %>%
  group_by(sample_code, healing_status) %>%
  summarise(
    fibroblast_cell_n = n(),
    across(all_of(score_cols), ~mean(.x, na.rm = TRUE)),
    mean_repair_axis = mean(repair_axis, na.rm = TRUE),
    mean_fibro_pseudotime = mean(fibro_pseudotime, na.rm = TRUE),
    .groups = "drop"
  )

score_tests <- bind_rows(lapply(c(score_cols, "mean_repair_axis", "mean_fibro_pseudotime"), function(feature) {
  dat <- dfu_sample_scores %>% select(healing_status, value = all_of(feature))
  data.frame(
    feature = feature,
    healer_n = sum(dat$healing_status == "Healer"),
    nonhealer_n = sum(dat$healing_status == "Non-healer"),
    healer_mean = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE) -
      mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(value ~ healing_status, data = dat)$p.value)
  )
})) %>%
  mutate(BH = p.adjust(p_wilcox, method = "BH"))

ref_variable_genes <- signature_rows %>%
  group_by(gene) %>%
  summarise(max_specificity = max(specificity), .groups = "drop") %>%
  arrange(desc(max_specificity)) %>%
  slice_head(n = 500) %>%
  pull(gene)
ref_variable_genes <- intersect(ref_variable_genes, rownames(dfu_log))
ref_avg_sub <- avg_log[ref_variable_genes, , drop = FALSE]

sample_avg <- sapply(unique(dfu_sample_scores$sample_code), function(sid) {
  cells <- rownames(fib_dfu@meta.data)[fib_dfu@meta.data$sample_code == sid & fib_dfu@meta.data$disease %in% c("DFU-healer", "DFU-nonhealer")]
  Matrix::rowMeans(dfu_log[ref_variable_genes, cells, drop = FALSE])
})

sample_condition_cor <- bind_rows(lapply(colnames(sample_avg), function(sid) {
  heal <- dfu_sample_scores$healing_status[dfu_sample_scores$sample_code == sid][1]
  data.frame(
    sample_code = sid,
    healing_status = heal,
    condition = colnames(ref_avg_sub),
    spearman_cor = apply(ref_avg_sub, 2, function(ref_vec) {
      suppressWarnings(cor(sample_avg[, sid], ref_vec, method = "spearman", use = "pairwise.complete.obs"))
    })
  )
}))

sample_condition_wide <- sample_condition_cor %>%
  tidyr::pivot_wider(names_from = condition, values_from = spearman_cor) %>%
  mutate(
    D7_minus_D0 = D7_wound - D0_skin,
    D30_minus_D0 = D30_wound - D0_skin,
    max_wound_minus_D0 = pmax(D1_wound, D7_wound, D30_wound, na.rm = TRUE) - D0_skin,
    best_reference_condition = conditions[max.col(as.matrix(select(., all_of(conditions))), ties.method = "first")]
  )

cor_tests <- bind_rows(lapply(c("D7_minus_D0", "D30_minus_D0", "max_wound_minus_D0", conditions), function(feature) {
  dat <- sample_condition_wide %>% select(healing_status, value = all_of(feature))
  data.frame(
    feature = feature,
    healer_mean = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE) -
      mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(value ~ healing_status, data = dat)$p.value)
  )
})) %>%
  mutate(BH = p.adjust(p_wilcox, method = "BH"))

write.table(fib_meta, file.path(out_dir, "GSE241132_fibroblast_metadata_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(signature_rows, file.path(out_dir, "GSE241132_fibroblast_condition_specific_signature_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wound_vs_skin_rows, file.path(out_dir, "GSE241132_fibroblast_wound_vs_skin_signature_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(avg_log) %>% mutate(gene = rownames(avg_log)), file.path(out_dir, "GSE241132_fibroblast_average_log_expression_by_condition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dfu_sample_scores, file.path(out_dir, "GSE165816_DFU_fibroblast_scores_against_GSE241132_wound_signatures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(score_tests, file.path(out_dir, "GSE165816_DFU_healer_vs_nonhealer_GSE241132_signature_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_condition_cor, file.path(out_dir, "GSE165816_DFU_fibroblast_reference_condition_correlations_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_condition_wide, file.path(out_dir, "GSE165816_DFU_fibroblast_reference_condition_correlations_wide.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cor_tests, file.path(out_dir, "GSE165816_DFU_fibroblast_reference_correlation_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_score_cols <- grep("D7_wound|D30_wound|D1_wound|D0_skin", score_cols, value = TRUE)
p1 <- dfu_sample_scores %>%
  select(sample_code, healing_status, fibroblast_cell_n, all_of(plot_score_cols)) %>%
  pivot_longer(all_of(plot_score_cols), names_to = "signature", values_to = "score") %>%
  ggplot(aes(x = healing_status, y = score, color = healing_status)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.12) +
  geom_point(aes(size = fibroblast_cell_n), position = position_jitter(width = 0.08, height = 0), alpha = 0.9) +
  facet_wrap(~signature, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Healer" = "#0072B2", "Non-healer" = "#D55E00")) +
  labs(x = NULL, y = "DFU fibroblast score", color = NULL, size = "Fibroblast cells") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "GSE165816_DFU_scores_against_GSE241132_wound_signatures.png"), p1, width = 13, height = 7, dpi = 240)

p2 <- sample_condition_cor %>%
  ggplot(aes(x = condition, y = sample_code, fill = spearman_cor)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_grid(healing_status ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#D55E00", mid = "grey95", high = "#0072B2", midpoint = 0) +
  labs(x = "GSE241132 reference condition", y = "GSE165816 DFU sample", fill = "Spearman\ncorrelation") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(fig_dir, "GSE165816_DFU_to_GSE241132_reference_condition_correlation_heatmap.png"), p2, width = 7, height = 6, dpi = 240)

p3 <- sample_condition_wide %>%
  select(sample_code, healing_status, D7_minus_D0, D30_minus_D0, max_wound_minus_D0) %>%
  pivot_longer(c(D7_minus_D0, D30_minus_D0, max_wound_minus_D0), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = healing_status, y = value, color = healing_status)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.12) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.2, alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Healer" = "#0072B2", "Non-healer" = "#D55E00")) +
  labs(x = NULL, y = "Reference wound alignment minus skin", color = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "GSE165816_DFU_reference_wound_alignment_metrics.png"), p3, width = 9, height = 4.8, dpi = 240)

message("Wrote human wound roadmap mapping outputs.")
