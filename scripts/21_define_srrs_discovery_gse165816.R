options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
obj_file <- file.path(project_dir, "results", "seurat_gse165816_foot_skin", "GSE165816_foot_skin_seurat_firstpass.rds")
d7_sig_file <- file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE241132_fibroblast_wound_vs_skin_signature_genes.tsv")
out_dir <- file.path(project_dir, "results", "srrs_framework")
fig_dir <- file.path(project_dir, "figures", "srrs_framework")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

module_sets <- list(
  fibroblast_repair_activation = c(
    "THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20",
    "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A",
    "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2"
  ),
  stromal_ligand_panel = c(
    "TNC", "IL11", "IL6", "INHBA", "THBS1", "SERPINE1", "CCL20",
    "WNT5A", "ADAM12", "PTGS2"
  ),
  vascular_perivascular_receiver_coupling = c(
    "ITGB1", "SDC4", "ITGA5", "ITGAV", "ITGA8", "CD47", "IL11RA",
    "IL6ST", "IL6R", "ENG", "ACVR1", "ACVR1B", "ACVR2A", "ACVR2B",
    "TGFBR3", "BAMBI", "LRP1", "PDGFRB", "KDR", "FLT1"
  ),
  negative_housekeeping_control = c(
    "ACTB", "GAPDH", "RPLP0", "RPS18", "B2M", "PPIA", "TBP", "HPRT1"
  )
)

if (!file.exists(d7_sig_file)) {
  stop("Missing D7 signature table. Run script 17_human_wound_roadmap_mapping_gse241132.R first.")
}
d7_tbl <- read.delim(d7_sig_file, check.names = FALSE)
d7_genes <- d7_tbl %>%
  filter(signature == "GSE241132_D7_wound_vs_D0_skin_fibroblast_up") %>%
  arrange(desc(delta_vs_skin)) %>%
  slice_head(n = 100) %>%
  pull(gene)
module_sets$d7_acute_wound_alignment <- d7_genes

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

z_by_discovery <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sx
}

wilcox_delta <- function(df, feature) {
  dat <- df %>% select(healing_status, value = all_of(feature)) %>% filter(is.finite(value))
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
}

exact_perm_p <- function(values, labels) {
  ok <- is.finite(values) & !is.na(labels)
  values <- values[ok]
  labels <- labels[ok]
  n_h <- sum(labels == "Healer")
  n <- length(values)
  obs <- mean(values[labels == "Healer"]) - mean(values[labels == "Non-healer"])
  cmb <- combn(seq_len(n), n_h)
  deltas <- apply(cmb, 2, function(idx) {
    mean(values[idx]) - mean(values[-idx])
  })
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

lodo_summary <- function(df, feature) {
  rows <- lapply(seq_len(nrow(df)), function(i) {
    sub <- df[-i, ]
    dat <- sub %>% select(healing_status, value = all_of(feature)) %>% filter(is.finite(value))
    data.frame(
      dropped_sample = df$sample_code[i],
      delta_healer_minus_nonhealer = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE) -
        mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE)
    )
  })
  bind_rows(rows)
}

message("Reading discovery Seurat object...")
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
receiver_cells <- meta %>%
  filter(broad_cell_type %in% c("endothelial", "pericyte_smc"), disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  pull(cell)

fib_expr <- expr[, fib_cells, drop = FALSE]
receiver_expr <- expr[, receiver_cells, drop = FALSE]
fib_meta <- meta[match(fib_cells, meta$cell), ]
receiver_meta <- meta[match(receiver_cells, meta$cell), ]

fib_scores <- data.frame(
  cell = fib_cells,
  sample_code = fib_meta$sample_code,
  fibroblast_repair_activation = score_cells(fib_expr, module_sets$fibroblast_repair_activation),
  stromal_ligand_panel = score_cells(fib_expr, module_sets$stromal_ligand_panel),
  d7_acute_wound_alignment = score_cells(fib_expr, module_sets$d7_acute_wound_alignment),
  negative_housekeeping_control = score_cells(fib_expr, module_sets$negative_housekeeping_control)
)

receiver_scores <- data.frame(
  cell = receiver_cells,
  sample_code = receiver_meta$sample_code,
  receiver_cell_type = receiver_meta$broad_cell_type,
  vascular_perivascular_receiver_coupling = score_cells(receiver_expr, module_sets$vascular_perivascular_receiver_coupling)
)

sample_components <- dfu_samples %>%
  left_join(
    fib_scores %>%
      group_by(sample_code) %>%
      summarise(
        fibroblast_cell_n = n(),
        fibroblast_repair_activation = mean(fibroblast_repair_activation, na.rm = TRUE),
        stromal_ligand_panel = mean(stromal_ligand_panel, na.rm = TRUE),
        d7_acute_wound_alignment = mean(d7_acute_wound_alignment, na.rm = TRUE),
        negative_housekeeping_control = mean(negative_housekeeping_control, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "sample_code"
  ) %>%
  left_join(
    receiver_scores %>%
      group_by(sample_code) %>%
      summarise(
        receiver_cell_n = n(),
        vascular_perivascular_receiver_coupling = mean(vascular_perivascular_receiver_coupling, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "sample_code"
  )

srrs_components <- c(
  "fibroblast_repair_activation",
  "stromal_ligand_panel",
  "vascular_perivascular_receiver_coupling",
  "d7_acute_wound_alignment"
)

standardisation <- bind_rows(lapply(srrs_components, function(feature) {
  data.frame(
    component = feature,
    discovery_mean = mean(sample_components[[feature]], na.rm = TRUE),
    discovery_sd = sd(sample_components[[feature]], na.rm = TRUE),
    weight = 0.25
  )
}))

for (feature in srrs_components) {
  sample_components[[paste0("z_", feature)]] <- z_by_discovery(sample_components[[feature]])
}
sample_components$SRRS <- rowMeans(sample_components[paste0("z_", srrs_components)], na.rm = TRUE)
sample_components$negative_housekeeping_control_z <- z_by_discovery(sample_components$negative_housekeeping_control)

test_features <- c(srrs_components, "SRRS", "negative_housekeeping_control_z")
tests <- bind_rows(lapply(test_features, function(feature) wilcox_delta(sample_components, feature))) %>%
  mutate(
    exact_permutation_p = sapply(feature, function(f) exact_perm_p(sample_components[[f]], sample_components$healing_status)),
    bootstrap_ci_low = sapply(feature, function(f) bootstrap_ci(sample_components[[f]], sample_components$healing_status)[1]),
    bootstrap_ci_high = sapply(feature, function(f) bootstrap_ci(sample_components[[f]], sample_components$healing_status)[2]),
    BH_wilcox = p.adjust(p_wilcox, method = "BH")
  )

lodo <- bind_rows(lapply(test_features, function(feature) {
  lodo_summary(sample_components, feature) %>% mutate(feature = feature)
})) %>%
  group_by(feature) %>%
  mutate(
    same_direction_as_full = sign(delta_healer_minus_nonhealer) ==
      sign(tests$delta_healer_minus_nonhealer[match(feature, tests$feature)])
  ) %>%
  ungroup()

lodo_summary_tbl <- lodo %>%
  group_by(feature) %>%
  summarise(
    lodo_min_delta = min(delta_healer_minus_nonhealer, na.rm = TRUE),
    lodo_max_delta = max(delta_healer_minus_nonhealer, na.rm = TRUE),
    lodo_same_direction = paste0(sum(same_direction_as_full, na.rm = TRUE), "/", n()),
    .groups = "drop"
  )
tests <- tests %>% left_join(lodo_summary_tbl, by = "feature")

gene_set_tbl <- bind_rows(lapply(names(module_sets), function(module) {
  data.frame(
    module = module,
    gene = module_sets[[module]],
    used_in_SRRS = module %in% srrs_components
  )
})) %>%
  group_by(module) %>%
  mutate(n_requested = n()) %>%
  ungroup()

write.table(gene_set_tbl, file.path(out_dir, "SRRS_locked_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(standardisation, file.path(out_dir, "SRRS_discovery_standardisation_parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_components, file.path(out_dir, "GSE165816_SRRS_discovery_sample_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tests, file.path(out_dir, "GSE165816_SRRS_discovery_healer_vs_nonhealer_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lodo, file.path(out_dir, "GSE165816_SRRS_discovery_leave_one_sample_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_components <- sample_components %>%
  select(sample_code, healing_status, all_of(test_features)) %>%
  pivot_longer(all_of(test_features), names_to = "feature", values_to = "value") %>%
  mutate(feature = factor(feature, levels = test_features)) %>%
  ggplot(aes(x = healing_status, y = value, color = healing_status)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.12) +
  geom_point(size = 2.4, position = position_jitter(width = 0.08, height = 0), alpha = 0.9) +
  facet_wrap(~feature, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Healer" = "#0072B2", "Non-healer" = "#D55E00")) +
  labs(x = NULL, y = "Sample-level score", color = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "GSE165816_SRRS_discovery_components.png"), p_components, width = 10.5, height = 7.2, dpi = 240)

p_lodo <- lodo %>%
  filter(feature %in% c("SRRS", srrs_components)) %>%
  ggplot(aes(x = delta_healer_minus_nonhealer, y = feature)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(color = "#0072B2", alpha = 0.8) +
  labs(x = "Leave-one-sample-out delta, healer minus non-healer", y = NULL) +
  theme_classic(base_size = 10)
ggsave(file.path(fig_dir, "GSE165816_SRRS_discovery_LODO.png"), p_lodo, width = 7, height = 4.5, dpi = 240)

message("Wrote SRRS discovery outputs.")
