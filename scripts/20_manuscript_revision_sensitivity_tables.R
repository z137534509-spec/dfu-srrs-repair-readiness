options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

set.seed(20260521)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_dir, "results", "manuscript_revision_burns_trauma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  read.delim(path, check.names = FALSE)
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

safe_wilcox <- function(value, group) {
  out <- tryCatch(
    suppressWarnings(wilcox.test(value ~ group)$p.value),
    error = function(e) NA_real_
  )
  out
}

delta_h_minus_nh <- function(value, group) {
  mean(value[group == "Healer"], na.rm = TRUE) - mean(value[group == "Non-healer"], na.rm = TRUE)
}

exact_permutation_p <- function(value, group) {
  group <- as.character(group)
  n <- length(group)
  n_healer <- sum(group == "Healer")
  obs <- delta_h_minus_nh(value, group)
  idx <- seq_len(n)
  combos <- combn(idx, n_healer)
  perm_delta <- apply(combos, 2, function(h_idx) {
    perm_group <- rep("Non-healer", n)
    perm_group[h_idx] <- "Healer"
    delta_h_minus_nh(value, perm_group)
  })
  p_two_sided <- mean(abs(perm_delta) >= abs(obs) - 1e-12)
  list(p = p_two_sided, n_permutations = length(perm_delta))
}

bootstrap_ci <- function(value, group, n_boot = 5000) {
  h <- value[group == "Healer"]
  n <- value[group == "Non-healer"]
  boots <- replicate(n_boot, mean(sample(h, length(h), replace = TRUE), na.rm = TRUE) -
    mean(sample(n, length(n), replace = TRUE), na.rm = TRUE))
  as.numeric(quantile(boots, c(0.025, 0.5, 0.975), na.rm = TRUE))
}

# Table 1: dataset and sample summaries --------------------------------------

obj <- readRDS(file.path(project_dir, "results", "seurat_gse165816_foot_skin", "GSE165816_foot_skin_seurat_firstpass.rds"))
meta <- obj@meta.data

sample_celltype <- meta %>%
  tibble::rownames_to_column("cell") %>%
  count(sample_code, geo_accession, disease, tissue, broad_cell_type, name = "cells") %>%
  pivot_wider(names_from = broad_cell_type, values_from = cells, values_fill = 0) %>%
  mutate(total_cells = rowSums(across(where(is.numeric))))

qc_by_sample <- read_tsv(file.path(project_dir, "results", "seurat_gse165816_foot_skin", "GSE165816_foot_skin_qc_by_sample.tsv"))
gse165816_sample_summary <- qc_by_sample %>%
  left_join(sample_celltype, by = c("sample_code", "geo_accession", "disease", "tissue")) %>%
  arrange(factor(disease, levels = c("DFU-healer", "DFU-nonhealer", "Non-DFU Diabetic", "Non-diabetic")), sample_code)
write.table(
  gse165816_sample_summary,
  file.path(out_dir, "Table1_GSE165816_sample_qc_celltype_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

gse165816_dataset <- meta %>%
  count(disease, name = "cells") %>%
  left_join(meta %>% distinct(sample_code, disease) %>% count(disease, name = "samples"), by = "disease") %>%
  mutate(summary = paste0(disease, ": ", samples, " samples, ", cells, " cells")) %>%
  summarise(
    dataset = "GSE165816",
    role = "Primary DFU outcome single-cell analysis",
    tissue_or_condition = "Human foot skin; DFU-healer, DFU-nonhealer, non-DFU diabetic and non-diabetic controls",
    samples = n_distinct(meta$sample_code),
    cells_or_profiles = nrow(meta),
    analysis_use = "Fibroblast repair-state analysis, pseudobulk DE, communication inference and in silico ligand prioritisation",
    sample_details = paste(summary, collapse = "; ")
  )

gse241132_meta <- read_tsv(file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE241132_fibroblast_metadata_used.tsv"))
gse241132_dataset <- gse241132_meta %>%
  count(timepoint, name = "cells") %>%
  left_join(gse241132_meta %>% distinct(orig.ident, timepoint) %>% count(timepoint, name = "samples"), by = "timepoint") %>%
  mutate(summary = paste0(timepoint, ": ", samples, " samples, ", cells, " fibroblasts")) %>%
  summarise(
    dataset = "GSE241132",
    role = "Physiological human acute wound fibroblast reference",
    tissue_or_condition = "Normal human skin and acute wounds at D0, D1, D7 and D30",
    samples = n_distinct(gse241132_meta$orig.ident),
    cells_or_profiles = nrow(gse241132_meta),
    analysis_use = "Reference alignment/contextualisation of DFU fibroblast states",
    sample_details = paste(summary, collapse = "; ")
  )

gse223964_file <- file.path(project_dir, "results", "validation_gse223964", "GSE223964_cell_counts_by_sample_cluster_celltype.tsv")
gse223964_dataset <- if (file.exists(gse223964_file)) {
  x <- read_tsv(gse223964_file)
  x %>%
    summarise(
      dataset = "GSE223964",
      role = "External diabetic vs non-diabetic scRNA trend analysis",
      tissue_or_condition = paste(sort(unique(condition)), collapse = "; "),
      samples = n_distinct(sample),
      cells_or_profiles = sum(cells, na.rm = TRUE),
      analysis_use = "Supportive cell-type and signature trend analysis; not an outcome validation cohort",
      sample_details = paste0("Samples: ", n_distinct(sample), "; annotated cells: ", sum(cells, na.rm = TRUE))
    )
} else {
  NULL
}

gse134431_file <- file.path(project_dir, "results", "bulk_validation_gse134431", "GSE134431_signature_scores.tsv")
gse134431_dataset <- if (file.exists(gse134431_file)) {
  x <- read_tsv(gse134431_file)
  x %>%
    summarise(
      dataset = "GSE134431",
      role = "Bulk transcriptomic supportive analysis",
      tissue_or_condition = paste(sort(unique(validation_group)), collapse = "; "),
      samples = n_distinct(sample_id),
      cells_or_profiles = n_distinct(sample_id),
      analysis_use = "Supportive signature-level analysis",
      sample_details = paste0("Validation groups: ", paste(sort(unique(validation_group)), collapse = "; "))
    )
} else {
  NULL
}

table1_dataset <- bind_rows(gse165816_dataset, gse241132_dataset, gse223964_dataset, gse134431_dataset)
write.table(table1_dataset, file.path(out_dir, "Table1_dataset_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Sensitivity analyses ---------------------------------------------------------

traj <- read_tsv(file.path(project_dir, "results", "fibroblast_trajectory_gse165816", "GSE165816_fibroblast_trajectory_sample_scores.tsv"))
wound_scores <- read_tsv(file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE165816_DFU_fibroblast_scores_against_GSE241132_wound_signatures.tsv"))
wound_corr <- read_tsv(file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE165816_DFU_fibroblast_reference_condition_correlations_wide.tsv"))

make_metric_long <- function(df, features, source) {
  df %>%
    select(sample_code, healing_status, all_of(features)) %>%
    pivot_longer(cols = all_of(features), names_to = "metric", values_to = "value") %>%
    mutate(source = source)
}

metric_long <- bind_rows(
  make_metric_long(
    traj,
    c("mean_repair_axis", "high_repair_fraction", "mean_he_fibro_score", "mean_ligand_score", "late_pseudotime_fraction"),
    "GSE165816 fibroblast repair-state"
  ),
  make_metric_long(
    wound_scores,
    c("GSE241132_D7_wound_vs_D0_skin_fibroblast_up", "GSE241132_D1_wound_vs_D0_skin_fibroblast_up"),
    "GSE241132 wound signature projection"
  ),
  make_metric_long(
    wound_corr,
    c("D7_wound", "D7_minus_D0", "max_wound_minus_D0", "D0_skin"),
    "GSE241132 reference correlation"
  )
)

metric_summary <- bind_rows(lapply(split(metric_long, metric_long$metric), function(dat) {
  dat <- dat %>% arrange(sample_code)
  perm <- exact_permutation_p(dat$value, dat$healing_status)
  ci <- bootstrap_ci(dat$value, dat$healing_status)
  data.frame(
    source = dat$source[1],
    metric = dat$metric[1],
    n_healer = sum(dat$healing_status == "Healer"),
    n_nonhealer = sum(dat$healing_status == "Non-healer"),
    healer_mean = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = delta_h_minus_nh(dat$value, dat$healing_status),
    bootstrap_delta_lcl95 = ci[1],
    bootstrap_delta_median = ci[2],
    bootstrap_delta_ucl95 = ci[3],
    p_wilcox = safe_wilcox(dat$value, dat$healing_status),
    p_exact_label_permutation = perm$p,
    n_exact_label_permutations = perm$n_permutations
  )
})) %>%
  mutate(BH_wilcox = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(p_exact_label_permutation, desc(abs(delta_healer_minus_nonhealer)))

write.table(
  metric_summary,
  file.path(out_dir, "GSE165816_core_metric_permutation_bootstrap_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

lodo_detail <- bind_rows(lapply(split(metric_long, metric_long$metric), function(dat) {
  obs_delta <- delta_h_minus_nh(dat$value, dat$healing_status)
  bind_rows(lapply(dat$sample_code, function(drop_sample) {
    keep <- dat %>% filter(sample_code != drop_sample)
    data.frame(
      source = dat$source[1],
      metric = dat$metric[1],
      dropped_sample = drop_sample,
      dropped_status = dat$healing_status[dat$sample_code == drop_sample],
      n_healer = sum(keep$healing_status == "Healer"),
      n_nonhealer = sum(keep$healing_status == "Non-healer"),
      delta_healer_minus_nonhealer = delta_h_minus_nh(keep$value, keep$healing_status),
      p_wilcox = safe_wilcox(keep$value, keep$healing_status),
      same_direction_as_observed = sign(delta_h_minus_nh(keep$value, keep$healing_status)) == sign(obs_delta)
    )
  }))
}))

lodo_summary <- lodo_detail %>%
  group_by(source, metric) %>%
  summarise(
    lodo_min_delta = min(delta_healer_minus_nonhealer, na.rm = TRUE),
    lodo_max_delta = max(delta_healer_minus_nonhealer, na.rm = TRUE),
    lodo_same_direction_n = sum(same_direction_as_observed, na.rm = TRUE),
    lodo_total_n = n(),
    lodo_all_same_direction = all(same_direction_as_observed, na.rm = TRUE),
    lodo_max_p_wilcox = max(p_wilcox, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(metric_summary %>% select(source, metric, observed_delta = delta_healer_minus_nonhealer, p_exact_label_permutation), by = c("source", "metric")) %>%
  arrange(p_exact_label_permutation)

write.table(lodo_detail, file.path(out_dir, "GSE165816_core_metric_leave_one_sample_out_detail.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lodo_summary, file.path(out_dir, "GSE165816_core_metric_leave_one_sample_out_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Consensus layer sensitivity --------------------------------------------------

lig <- read_tsv(file.path(project_dir, "results", "fibroblast_communication_consensus_gse165816", "GSE165816_fibroblast_communication_ligand_consensus.tsv"))
axis <- read_tsv(file.path(project_dir, "results", "fibroblast_communication_consensus_gse165816", "GSE165816_fibroblast_communication_axis_consensus.tsv"))

component_weights <- c(
  fib_component = 1.4,
  integrated_component = 1.2,
  nichenet_component = 1.2,
  lr_component = 1.2,
  virtual_component = 1.4,
  repair_axis_component = 1.0,
  normal_wound_component = 1.0,
  tf_component = 0.8
)

score_without <- function(df, drop_component = NA_character_) {
  weights <- component_weights
  if (!is.na(drop_component)) weights <- weights[names(weights) != drop_component]
  raw <- as.matrix(df[, names(weights), drop = FALSE]) %*% weights
  score <- 100 * as.numeric(raw) / sum(weights)
  rank <- rank(-score, ties.method = "min")
  data.frame(gene = df$gene, drop_component = ifelse(is.na(drop_component), "none", drop_component), sensitivity_score = score, sensitivity_rank = rank)
}

consensus_sensitivity_long <- bind_rows(
  score_without(lig, NA_character_),
  bind_rows(lapply(names(component_weights), function(comp) score_without(lig, comp)))
)

consensus_sensitivity_summary <- consensus_sensitivity_long %>%
  group_by(gene) %>%
  summarise(
    full_rank = sensitivity_rank[drop_component == "none"],
    full_score = sensitivity_score[drop_component == "none"],
    min_rank_after_dropping_one_layer = min(sensitivity_rank[drop_component != "none"], na.rm = TRUE),
    max_rank_after_dropping_one_layer = max(sensitivity_rank[drop_component != "none"], na.rm = TRUE),
    median_rank_after_dropping_one_layer = median(sensitivity_rank[drop_component != "none"], na.rm = TRUE),
    always_top5_after_dropping_one_layer = all(sensitivity_rank[drop_component != "none"] <= 5, na.rm = TRUE),
    always_top10_after_dropping_one_layer = all(sensitivity_rank[drop_component != "none"] <= 10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lig %>% select(gene, consensus_tier, evidence_count, consensus_score), by = "gene") %>%
  arrange(full_rank)

write.table(consensus_sensitivity_long, file.path(out_dir, "Ligand_consensus_leave_one_evidence_layer_out_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(consensus_sensitivity_summary, file.path(out_dir, "Ligand_consensus_leave_one_evidence_layer_out_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Table 2: candidate ligand axis summary ---------------------------------------

top_ligands <- c("TNC", "IL11", "IL6", "INHBA", "THBS1", "SERPINE1", "CCL20")

interpretation <- c(
  TNC = "Matrix organisation and fibroblast-endothelial/perivascular cell-matrix communication",
  IL11 = "gp130-related fibroblast repair signalling",
  IL6 = "gp130-related inflammatory repair signalling",
  INHBA = "TGF-like stromal and endothelial/perivascular coupling",
  THBS1 = "Matrix remodelling, integrin/CD47 signalling and protease balance",
  SERPINE1 = "Protease balance and matrix remodelling",
  CCL20 = "Inflammatory-repair chemokine signalling"
)

validation_priority <- c(
  TNC = "High",
  IL11 = "High",
  IL6 = "High",
  INHBA = "High",
  THBS1 = "Medium-high",
  SERPINE1 = "Medium-high",
  CCL20 = "Medium"
)

table2 <- lig %>%
  filter(gene %in% top_ligands) %>%
  left_join(
    axis %>%
      filter(ligand %in% top_ligands) %>%
      group_by(ligand) %>%
      arrange(desc(axis_consensus_score), .by_group = TRUE) %>%
      summarise(
        top_lr_axes = collapse_unique(head(interaction, 4)),
        top_receiver_compartments = collapse_unique(head(receiver, 4)),
        best_axis_score = first(axis_consensus_score),
        best_axis_p_wilcox = first(p_wilcox),
        .groups = "drop"
      ),
    by = c("gene" = "ligand")
  ) %>%
  mutate(
    fibroblast_DE = paste0("log2FC ", round(fib_log2fc, 2), ", BH ", signif(fib_padj, 2)),
    repair_axis = paste0("rho ", round(repair_axis_spearman_rho, 2)),
    acute_wound_context = ifelse(is.na(wound_conditions_up), "Not in D1/D7/D30 wound-up list", wound_conditions_up),
    regulon_context = ifelse(is.na(tf_leading_edge_pathways), "Not leading-edge supported", tf_leading_edge_pathways),
    clinical_interpretation = unname(interpretation[gene]),
    proposed_validation_priority = unname(validation_priority[gene])
  ) %>%
  select(
    ligand = gene,
    consensus_tier,
    evidence_count,
    consensus_score,
    fibroblast_DE,
    nichenet_max_pearson,
    top_lr_axes,
    top_receiver_compartments,
    virtual_best_receiver,
    virtual_best_priority,
    repair_axis,
    acute_wound_context,
    regulon_context,
    clinical_interpretation,
    proposed_validation_priority
  ) %>%
  arrange(desc(consensus_score))

write.table(table2, file.path(out_dir, "Table2_candidate_ligand_axis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Wrote manuscript revision outputs to: ", out_dir)
