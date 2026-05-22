options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
fib_file <- file.path(project_dir, "results", "fibroblast_trajectory_gse165816", "GSE165816_fibroblast_trajectory_seurat.rds")
dorothea_file <- file.path(project_dir, "data", "metadata", "omnipath_dorothea_interactions_2026-05-21.tsv")
out_dir <- file.path(project_dir, "results", "fibroblast_regulon_dorothea_gse165816")
fig_dir <- file.path(project_dir, "figures", "fibroblast_regulon_dorothea_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fib <- readRDS(fib_file)
DefaultAssay(fib) <- "RNA"
expr <- LayerData(fib, assay = "RNA", layer = "data")
counts <- LayerData(fib, assay = "RNA", layer = "counts")

meta <- fib@meta.data %>%
  mutate(cell = rownames(.))

q25 <- quantile(meta$repair_axis, 0.25, na.rm = TRUE)
q75 <- quantile(meta$repair_axis, 0.75, na.rm = TRUE)
meta$repair_state_group <- case_when(
  meta$repair_axis <= q25 ~ "low_repair",
  meta$repair_axis >= q75 ~ "high_repair",
  TRUE ~ "middle"
)
fib@meta.data$repair_state_group <- meta$repair_state_group

high_cells <- meta$cell[meta$repair_state_group == "high_repair"]
low_cells <- meta$cell[meta$repair_state_group == "low_repair"]

avg_high <- Matrix::rowMeans(expr[, high_cells, drop = FALSE])
avg_low <- Matrix::rowMeans(expr[, low_cells, drop = FALSE])
pct_high <- Matrix::rowMeans(counts[, high_cells, drop = FALSE] > 0)
pct_low <- Matrix::rowMeans(counts[, low_cells, drop = FALSE] > 0)

repair_delta <- data.frame(
  gene = rownames(expr),
  avg_log_high = avg_high,
  avg_log_low = avg_low,
  delta_high_minus_low = avg_high - avg_low,
  pct_high = pct_high,
  pct_low = pct_low
) %>%
  mutate(max_pct = pmax(pct_high, pct_low)) %>%
  filter(max_pct >= 0.02) %>%
  arrange(desc(delta_high_minus_low))

rank_stats <- repair_delta$delta_high_minus_low
names(rank_stats) <- repair_delta$gene
rank_stats <- sort(rank_stats, decreasing = TRUE)
rank_stats <- rank_stats[is.finite(rank_stats)]

dor <- read.delim(dorothea_file, check.names = FALSE)
dor_use <- dor %>%
  filter(
    dorothea_level %in% c("A", "B", "C"),
    source_genesymbol != "",
    target_genesymbol != "",
    target_genesymbol %in% names(rank_stats)
  ) %>%
  distinct(source_genesymbol, target_genesymbol, dorothea_level, consensus_stimulation, consensus_inhibition, .keep_all = TRUE)

regulons <- split(dor_use$target_genesymbol, dor_use$source_genesymbol)
regulons <- lapply(regulons, unique)
regulons <- regulons[lengths(regulons) >= 10 & lengths(regulons) <= 1500]

fg <- fgsea(
  pathways = regulons,
  stats = rank_stats,
  minSize = 10,
  maxSize = 1500,
  eps = 0
) %>%
  as.data.frame() %>%
  arrange(padj, desc(NES)) %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1)))

candidate_tfs <- c(
  "STAT3", "HIF1A", "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
  "CEBPB", "CEBPD", "SMAD2", "SMAD3", "SMAD4", "NFKB1", "RELA", "REL",
  "TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "EGR1", "KLF6", "KLF4", "ATF3"
)
candidate_tfs <- intersect(candidate_tfs, names(regulons))

score_gene_set <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 5) return(rep(NA_real_, ncol(expr_mat)))
  mat <- as.matrix(expr_mat[genes, , drop = FALSE])
  keep <- apply(mat, 1, sd, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 5) return(rep(NA_real_, ncol(expr_mat)))
  z <- t(scale(t(mat)))
  colMeans(z, na.rm = TRUE)
}

tf_scores <- as.data.frame(sapply(candidate_tfs, function(tf) score_gene_set(expr, regulons[[tf]])))
tf_scores$cell <- rownames(tf_scores)
tf_scores <- tf_scores %>%
  left_join(meta %>% select(cell, sample_code, disease, repair_axis, repair_state_group), by = "cell")

tf_sample_scores <- tf_scores %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  mutate(healing_status = ifelse(disease == "DFU-healer", "Healer", "Non-healer")) %>%
  group_by(sample_code, healing_status) %>%
  summarise(
    fibroblast_cell_n = n(),
    across(all_of(candidate_tfs), ~mean(.x, na.rm = TRUE)),
    mean_repair_axis = mean(repair_axis, na.rm = TRUE),
    .groups = "drop"
  )

tf_sample_tests <- bind_rows(lapply(candidate_tfs, function(tf) {
  dat <- tf_sample_scores %>% select(healing_status, value = all_of(tf))
  data.frame(
    TF = tf,
    target_n = length(regulons[[tf]]),
    healer_mean = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE) -
      mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(value ~ healing_status, data = dat)$p.value)
  )
})) %>%
  left_join(fg %>% select(pathway, NES, padj, leadingEdge), by = c("TF" = "pathway")) %>%
  mutate(BH_sample = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(desc(delta_healer_minus_nonhealer))

write.table(repair_delta, file.path(out_dir, "GSE165816_fibroblast_high_vs_low_repair_delta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dor_use, file.path(out_dir, "DoRothEA_ABC_regulon_edges_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fg, file.path(out_dir, "GSE165816_fibroblast_high_repair_DoRothEA_fgsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tf_sample_scores, file.path(out_dir, "GSE165816_fibroblast_candidate_TF_regulon_scores_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tf_sample_tests, file.path(out_dir, "GSE165816_fibroblast_candidate_TF_regulon_healer_vs_nonhealer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_fg <- fg %>%
  filter(NES > 0, padj < 0.05) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

if (nrow(top_fg) > 0) {
  p1 <- ggplot(top_fg, aes(x = pathway, y = NES, fill = padj)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_gradient(low = "#0072B2", high = "#D55E00", trans = "reverse") +
    labs(x = NULL, y = "FGSEA NES in high-repair fibroblasts", fill = "padj") +
    theme_classic(base_size = 10)
  ggsave(file.path(fig_dir, "GSE165816_fibroblast_high_repair_top_DoRothEA_TF_fgsea.png"), p1, width = 7.5, height = 7, dpi = 240)
}

plot_tfs <- tf_sample_tests %>%
  filter(!is.na(delta_healer_minus_nonhealer)) %>%
  arrange(desc(delta_healer_minus_nonhealer)) %>%
  slice_head(n = 16) %>%
  pull(TF)

if (length(plot_tfs) > 0) {
  p2 <- tf_sample_scores %>%
    select(sample_code, healing_status, fibroblast_cell_n, all_of(plot_tfs)) %>%
    pivot_longer(all_of(plot_tfs), names_to = "TF", values_to = "regulon_score") %>%
    ggplot(aes(x = healing_status, y = regulon_score, color = healing_status)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.12) +
    geom_point(aes(size = fibroblast_cell_n), position = position_jitter(width = 0.08, height = 0), alpha = 0.9) +
    facet_wrap(~TF, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c("Healer" = "#0072B2", "Non-healer" = "#D55E00")) +
    labs(x = NULL, y = "Candidate TF regulon score", color = NULL, size = "Fibroblast cells") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(fig_dir, "GSE165816_fibroblast_candidate_TF_regulon_scores_by_sample.png"), p2, width = 10, height = 7, dpi = 240)
}

message("Wrote fibroblast regulon outputs.")
