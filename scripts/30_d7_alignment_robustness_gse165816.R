options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
obj_file <- file.path(project_dir, "results", "seurat_gse165816_foot_skin", "GSE165816_foot_skin_seurat_firstpass.rds")
sig_file <- file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE241132_fibroblast_wound_vs_skin_signature_genes.tsv")
avg_file <- file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE241132_fibroblast_average_log_expression_by_condition.tsv")
out_dir <- file.path(project_dir, "results", "srrs_framework", "d7_alignment_robustness")
fig_dir <- file.path(project_dir, "figures", "srrs_framework", "d7_alignment_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

score_cells <- function(expr, genes) {
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
  zmat <- t(scale(t(mat)))
  colMeans(zmat, na.rm = TRUE)
}

exact_perm_p <- function(values, labels) {
  ok <- is.finite(values) & !is.na(labels)
  values <- values[ok]
  labels <- labels[ok]
  n_h <- sum(labels == "Healer")
  n <- length(values)
  obs <- mean(values[labels == "Healer"]) - mean(values[labels == "Non-healer"])
  cmb <- combn(seq_len(n), n_h)
  deltas <- apply(cmb, 2, function(idx) mean(values[idx]) - mean(values[-idx]))
  (sum(abs(deltas) >= abs(obs)) + 1) / (length(deltas) + 1)
}

bootstrap_ci <- function(values, labels, n_boot = 5000, seed = 1) {
  set.seed(seed)
  h <- values[labels == "Healer" & is.finite(values)]
  n <- values[labels == "Non-healer" & is.finite(values)]
  if (length(h) < 2 || length(n) < 2) return(c(NA_real_, NA_real_))
  deltas <- replicate(n_boot, mean(sample(h, replace = TRUE)) - mean(sample(n, replace = TRUE)))
  as.numeric(quantile(deltas, c(0.025, 0.975), na.rm = TRUE))
}

test_feature <- function(df, feature) {
  dat <- df %>% select(healing_status, value = all_of(feature)) %>% filter(is.finite(value))
  ci <- bootstrap_ci(dat$value, dat$healing_status)
  data.frame(
    feature = feature,
    healer_n = sum(dat$healing_status == "Healer"),
    nonhealer_n = sum(dat$healing_status == "Non-healer"),
    healer_mean = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE) -
      mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(value ~ healing_status, data = dat)$p.value),
    exact_permutation_p = exact_perm_p(dat$value, dat$healing_status),
    bootstrap_ci_low = ci[1],
    bootstrap_ci_high = ci[2]
  )
}

lodo_summary <- function(df, feature) {
  full_delta <- test_feature(df, feature)$delta_healer_minus_nonhealer
  rows <- lapply(seq_len(nrow(df)), function(i) {
    sub <- df[-i, ]
    delta <- test_feature(sub, feature)$delta_healer_minus_nonhealer
    data.frame(
      feature = feature,
      dropped_sample = df$sample_code[i],
      delta_healer_minus_nonhealer = delta,
      same_direction_as_full = sign(delta) == sign(full_delta)
    )
  })
  bind_rows(rows)
}

top_genes <- function(tbl, score_col, n = 100) {
  tbl %>%
    filter(is.finite(.data[[score_col]]), .data[[score_col]] > 0) %>%
    arrange(desc(.data[[score_col]])) %>%
    slice_head(n = n) %>%
    pull(gene)
}

message("Building alternative acute wound alignment gene sets...")
sig <- read.delim(sig_file, check.names = FALSE)
avg <- read.delim(avg_file, check.names = FALSE) %>%
  mutate(
    delta_D1_D0 = D1_wound - D0_skin,
    delta_D7_D0 = D7_wound - D0_skin,
    delta_D30_D0 = D30_wound - D0_skin,
    delta_D7_D1 = D7_wound - D1_wound,
    delta_D7_D30 = D7_wound - D30_wound,
    delta_D7D30_D0 = ((D7_wound + D30_wound) / 2) - D0_skin,
    D7_consensus_score = pmin(delta_D7_D0, delta_D7_D1, delta_D7_D30)
  )

gene_sets <- list(
  D1_vs_D0_top100 = sig %>%
    filter(signature == "GSE241132_D1_wound_vs_D0_skin_fibroblast_up") %>%
    arrange(desc(delta_vs_skin)) %>%
    slice_head(n = 100) %>%
    pull(gene),
  D7_vs_D0_top100_original = sig %>%
    filter(signature == "GSE241132_D7_wound_vs_D0_skin_fibroblast_up") %>%
    arrange(desc(delta_vs_skin)) %>%
    slice_head(n = 100) %>%
    pull(gene),
  D30_vs_D0_top100 = sig %>%
    filter(signature == "GSE241132_D30_wound_vs_D0_skin_fibroblast_up") %>%
    arrange(desc(delta_vs_skin)) %>%
    slice_head(n = 100) %>%
    pull(gene),
  D7_vs_D1_top100 = top_genes(avg, "delta_D7_D1", 100),
  D7_vs_D30_top100 = top_genes(avg, "delta_D7_D30", 100),
  D7D30_vs_D0_top100 = top_genes(avg, "delta_D7D30_D0", 100),
  D7_consensus_vs_D0_D1_D30_top100 = top_genes(avg, "D7_consensus_score", 100)
)

gene_set_tbl <- bind_rows(lapply(names(gene_sets), function(name) {
  data.frame(signature = name, gene = gene_sets[[name]], n_genes = length(gene_sets[[name]]))
}))
write.table(gene_set_tbl, file.path(out_dir, "GSE241132_alternative_alignment_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Reading GSE165816 discovery object...")
obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"
expr <- LayerData(obj, assay = "RNA", layer = "data")
meta <- obj@meta.data %>% mutate(cell = rownames(.))

dfu_samples <- meta %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  distinct(sample_code, disease) %>%
  mutate(healing_status = ifelse(disease == "DFU-healer", "Healer", "Non-healer")) %>%
  arrange(sample_code)

fib_cells <- meta %>%
  filter(broad_cell_type == "fibroblast_stromal", disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  pull(cell)
fib_expr <- expr[, fib_cells, drop = FALSE]
fib_meta <- meta[match(fib_cells, meta$cell), ]

cell_scores <- data.frame(cell = fib_cells, sample_code = fib_meta$sample_code)
for (name in names(gene_sets)) {
  cell_scores[[name]] <- score_cells(fib_expr, gene_sets[[name]])
}

sample_scores <- dfu_samples %>%
  left_join(
    cell_scores %>%
      group_by(sample_code) %>%
      summarise(across(all_of(names(gene_sets)), ~ mean(.x, na.rm = TRUE)), .groups = "drop"),
    by = "sample_code"
  )

tests <- bind_rows(lapply(names(gene_sets), function(feature) test_feature(sample_scores, feature))) %>%
  mutate(BH_wilcox = p.adjust(p_wilcox, method = "BH"))
lodo <- bind_rows(lapply(names(gene_sets), function(feature) lodo_summary(sample_scores, feature)))
lodo_sum <- lodo %>%
  group_by(feature) %>%
  summarise(
    lodo_min_delta = min(delta_healer_minus_nonhealer, na.rm = TRUE),
    lodo_max_delta = max(delta_healer_minus_nonhealer, na.rm = TRUE),
    lodo_same_direction = paste0(sum(same_direction_as_full, na.rm = TRUE), "/", n()),
    .groups = "drop"
  )
tests <- tests %>% left_join(lodo_sum, by = "feature")

write.table(sample_scores, file.path(out_dir, "GSE165816_alternative_alignment_sample_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tests, file.path(out_dir, "GSE165816_alternative_alignment_healer_nonhealer_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lodo, file.path(out_dir, "GSE165816_alternative_alignment_leave_one_sample_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- tests %>%
  arrange(delta_healer_minus_nonhealer) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  ggplot(aes(x = delta_healer_minus_nonhealer, y = feature)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbarh(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), height = 0.15, color = "grey45") +
  geom_point(color = "#0072B2", size = 2.4) +
  labs(x = "Healer minus non-healer delta, bootstrap 95% CI", y = NULL) +
  theme_classic(base_size = 10)
ggsave(file.path(fig_dir, "GSE165816_D7_alignment_robustness_effects.png"), p, width = 7.5, height = 4.8, dpi = 240)

message("Wrote D7 alignment robustness outputs.")
