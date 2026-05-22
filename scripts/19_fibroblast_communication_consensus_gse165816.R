options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_dir, "results", "fibroblast_communication_consensus_gse165816")
fig_dir <- file.path(project_dir, "figures", "fibroblast_communication_consensus_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  read.delim(path, check.names = FALSE)
}

rescale01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  rng <- range(x[ok], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(diff(rng)) < 1e-12) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (x[ok] - rng[1]) / diff(rng)
  }
  out
}

split_components <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  unique(unlist(strsplit(x, "[_;|;/, ]+")))
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

integrated_file <- file.path(project_dir, "results", "integrated_candidates_gse165816", "GSE165816_integrated_fibroblast_sender_ligand_priority.tsv")
shortlist_file <- file.path(project_dir, "results", "integrated_candidates_gse165816", "GSE165816_integrated_fibroblast_sender_ligand_shortlist.tsv")
nichenet_file <- file.path(project_dir, "results", "nichenet_gse165816", "GSE165816_nichenet_ligand_activities.tsv")
lr_axis_file <- file.path(project_dir, "results", "omnipath_lr_gse165816", "GSE165816_omnipath_lr_candidate_axis_tests.tsv")
virtual_file <- file.path(project_dir, "results", "virtual_perturbation_gse165816", "GSE165816_virtual_ligand_knockdown_summary.tsv")
all_de_file <- file.path(project_dir, "results", "pseudobulk_de_gse165816", "GSE165816_all_focus_celltypes_DESeq2_healer_vs_nonhealer.tsv")
trajectory_file <- file.path(project_dir, "results", "fibroblast_trajectory_gse165816", "GSE165816_fibroblast_candidate_gene_pseudotime_correlations.tsv")
wound_sig_file <- file.path(project_dir, "results", "human_wound_roadmap_gse241132", "GSE241132_fibroblast_wound_vs_skin_signature_genes.tsv")
tf_fgsea_file <- file.path(project_dir, "results", "fibroblast_regulon_dorothea_gse165816", "GSE165816_fibroblast_high_repair_DoRothEA_fgsea.tsv")

integrated <- read_tsv(integrated_file)
shortlist <- read_tsv(shortlist_file)
nichenet <- read_tsv(nichenet_file)
lr_axis <- read_tsv(lr_axis_file)
virtual <- read_tsv(virtual_file)
all_de <- read_tsv(all_de_file)
trajectory <- read_tsv(trajectory_file)
wound_sig <- read_tsv(wound_sig_file)
tf_fgsea <- read_tsv(tf_fgsea_file)

focus_receivers <- c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte", "pericyte_smc")

candidate_ligands <- sort(unique(c(
  shortlist$gene,
  integrated$gene[integrated$priority_score >= 6],
  virtual$ligand
)))
candidate_ligands <- candidate_ligands[!is.na(candidate_ligands) & nzchar(candidate_ligands)]

fib_de <- all_de %>%
  filter(cell_type == "fibroblast_stromal") %>%
  select(gene, fib_log2fc = log2FoldChange, fib_padj = padj, fib_baseMean = baseMean)

niche_summary <- nichenet %>%
  filter(sender == "fibroblast_stromal", receiver %in% focus_receivers) %>%
  group_by(test_ligand) %>%
  arrange(desc(pearson), .by_group = TRUE) %>%
  summarise(
    nichenet_max_pearson = max(pearson, na.rm = TRUE),
    nichenet_mean_pearson = mean(pearson, na.rm = TRUE),
    nichenet_best_receiver = first(receiver),
    nichenet_best_auroc = first(auroc),
    nichenet_positive_receiver_n = sum(pearson > 0, na.rm = TRUE),
    nichenet_positive_receivers = collapse_unique(receiver[pearson > 0]),
    .groups = "drop"
  ) %>%
  rename(gene = test_ligand)

lr_summary <- lr_axis %>%
  filter(sender == "fibroblast_stromal", receiver %in% focus_receivers, delta_healer_minus_nonhealer > 0) %>%
  group_by(ligand) %>%
  arrange(p_wilcox, desc(ratio_healer_over_nonhealer), .by_group = TRUE) %>%
  summarise(
    lr_best_interaction = first(interaction),
    lr_best_receiver = first(receiver),
    lr_best_raw_p = first(p_wilcox),
    lr_best_candidate_padj = first(padj_candidate_bh),
    lr_best_ratio = first(ratio_healer_over_nonhealer),
    lr_healer_higher_pairs_raw_p05 = sum(p_wilcox < 0.05, na.rm = TRUE),
    lr_healer_higher_pairs = n(),
    lr_receivers = collapse_unique(receiver),
    .groups = "drop"
  ) %>%
  rename(gene = ligand)

virtual_summary <- virtual %>%
  group_by(ligand) %>%
  arrange(desc(perturbation_priority), .by_group = TRUE) %>%
  summarise(
    virtual_best_receiver = first(receiver),
    virtual_best_priority = first(perturbation_priority),
    virtual_best_knockdown_shift = first(predicted_knockdown_shift_toward_nonhealer),
    virtual_receiver_n_priority10 = sum(perturbation_priority >= 10, na.rm = TRUE),
    virtual_receivers_priority10 = collapse_unique(receiver[perturbation_priority >= 10]),
    .groups = "drop"
  ) %>%
  rename(gene = ligand)

wound_summary <- wound_sig %>%
  filter(gene %in% candidate_ligands) %>%
  group_by(gene) %>%
  summarise(
    wound_d1_up = any(condition == "D1_wound"),
    wound_d7_up = any(condition == "D7_wound"),
    wound_d30_up = any(condition == "D30_wound"),
    wound_max_delta_vs_skin = max(delta_vs_skin, na.rm = TRUE),
    wound_conditions_up = collapse_unique(condition),
    .groups = "drop"
  )

top_tf_pathways <- c("RELA", "NFKB1", "NFKB2", "REL", "STAT3", "SMAD3", "SMAD4", "HIF1A", "EPAS1", "JUN", "FOS", "CEBPB", "EGR1")
tf_edges <- do.call(rbind, lapply(seq_len(nrow(tf_fgsea)), function(i) {
  genes <- split_components(tf_fgsea$leadingEdge[i])
  if (length(genes) == 0) return(NULL)
  data.frame(
    pathway = tf_fgsea$pathway[i],
    gene = genes,
    NES = tf_fgsea$NES[i],
    padj = tf_fgsea$padj[i]
  )
}))

tf_summary <- tf_edges %>%
  filter(pathway %in% top_tf_pathways, NES > 0, padj < 0.05, gene %in% candidate_ligands) %>%
  group_by(gene) %>%
  summarise(
    tf_leading_edge_n = n_distinct(pathway),
    tf_leading_edge_pathways = collapse_unique(pathway),
    .groups = "drop"
  )

ligand_consensus <- data.frame(gene = candidate_ligands) %>%
  left_join(integrated %>% select(gene, integrated_priority_score = priority_score), by = "gene") %>%
  left_join(fib_de, by = "gene") %>%
  left_join(niche_summary, by = "gene") %>%
  left_join(lr_summary, by = "gene") %>%
  left_join(virtual_summary, by = "gene") %>%
  left_join(trajectory %>% select(gene, repair_axis_spearman_rho = spearman_rho, repair_axis_BH = BH), by = "gene") %>%
  left_join(wound_summary, by = "gene") %>%
  left_join(tf_summary, by = "gene")

ligand_consensus <- ligand_consensus %>%
  mutate(
    wound_d1_up = ifelse(is.na(wound_d1_up), FALSE, wound_d1_up),
    wound_d7_up = ifelse(is.na(wound_d7_up), FALSE, wound_d7_up),
    wound_d30_up = ifelse(is.na(wound_d30_up), FALSE, wound_d30_up),
    tf_leading_edge_n = ifelse(is.na(tf_leading_edge_n), 0, tf_leading_edge_n),
    fib_de_support = !is.na(fib_log2fc) & fib_log2fc > 0 & !is.na(fib_padj) & fib_padj < 0.1,
    integrated_support = !is.na(integrated_priority_score) & integrated_priority_score >= 6,
    nichenet_support = !is.na(nichenet_max_pearson) & nichenet_max_pearson > 0.04,
    omnipath_lr_support = !is.na(lr_best_raw_p) & lr_best_raw_p < 0.05,
    virtual_perturbation_support = !is.na(virtual_best_priority) & virtual_best_priority >= 10,
    repair_axis_support = !is.na(repair_axis_spearman_rho) & repair_axis_spearman_rho > 0.15 & !is.na(repair_axis_BH) & repair_axis_BH < 0.1,
    normal_wound_support = wound_d7_up | wound_d1_up,
    tf_regulon_support = tf_leading_edge_n >= 2,
    evidence_count = fib_de_support + integrated_support + nichenet_support + omnipath_lr_support +
      virtual_perturbation_support + repair_axis_support + normal_wound_support + tf_regulon_support
  )

ligand_consensus <- ligand_consensus %>%
  mutate(
    fib_component = ifelse(fib_log2fc > 0, rescale01(pmin(-log10(pmax(fib_padj, 1e-300)), 10) * pmax(fib_log2fc, 0)), 0),
    integrated_component = rescale01(integrated_priority_score),
    nichenet_component = rescale01(nichenet_max_pearson),
    lr_component = rescale01(pmin(-log10(pmax(lr_best_raw_p, 1e-300)), 10)),
    virtual_component = rescale01(virtual_best_priority),
    repair_axis_component = rescale01(pmax(repair_axis_spearman_rho, 0)),
    normal_wound_component = pmin(as.numeric(wound_d1_up) + as.numeric(wound_d7_up) + 0.5 * as.numeric(wound_d30_up), 2) / 2,
    tf_component = pmin(tf_leading_edge_n, 4) / 4,
    consensus_score = 100 * (
      1.4 * fib_component +
        1.2 * integrated_component +
        1.2 * nichenet_component +
        1.2 * lr_component +
        1.4 * virtual_component +
        1.0 * repair_axis_component +
        1.0 * normal_wound_component +
        0.8 * tf_component
    ) / 9.2,
    consensus_tier = case_when(
      evidence_count >= 7 ~ "Tier 1: convergent",
      evidence_count >= 6 ~ "Tier 2: strong",
      evidence_count >= 5 ~ "Tier 3: supportive",
      TRUE ~ "Exploratory"
    )
  ) %>%
  arrange(desc(consensus_score), desc(evidence_count), desc(virtual_best_priority))

ligand_out <- file.path(out_dir, "GSE165816_fibroblast_communication_ligand_consensus.tsv")
write.table(ligand_consensus, ligand_out, sep = "\t", quote = FALSE, row.names = FALSE)

evidence_layers <- c(
  fib_de_support = "Fibroblast DE",
  integrated_support = "Prior shortlist",
  nichenet_support = "NicheNet",
  omnipath_lr_support = "OmniPath LR",
  virtual_perturbation_support = "Virtual perturbation",
  repair_axis_support = "Repair axis",
  normal_wound_support = "Normal wound",
  tf_regulon_support = "TF regulon"
)

evidence_long <- do.call(rbind, lapply(names(evidence_layers), function(col) {
  data.frame(
    gene = ligand_consensus$gene,
    consensus_score = ligand_consensus$consensus_score,
    consensus_tier = ligand_consensus$consensus_tier,
    layer = evidence_layers[[col]],
    supported = as.logical(ligand_consensus[[col]])
  )
}))
evidence_long$layer <- factor(evidence_long$layer, levels = unname(evidence_layers))
evidence_long$gene <- factor(evidence_long$gene, levels = rev(ligand_consensus$gene))

evidence_out <- file.path(out_dir, "GSE165816_fibroblast_communication_ligand_evidence_matrix_long.tsv")
write.table(evidence_long, evidence_out, sep = "\t", quote = FALSE, row.names = FALSE)

p1 <- ggplot(evidence_long, aes(x = layer, y = gene, fill = supported)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(`FALSE` = "#E7ECEF", `TRUE` = "#2E6F95")) +
  labs(x = NULL, y = NULL, fill = "Supported", title = "Fibroblast-sender ligand consensus evidence") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(fig_dir, "GSE165816_fibroblast_ligand_consensus_evidence_heatmap.png"), p1, width = 8.5, height = max(4, 0.28 * nrow(ligand_consensus) + 1.8), dpi = 300)

receptor_de_lookup <- all_de %>%
  filter(cell_type %in% focus_receivers) %>%
  select(receiver = cell_type, receptor_gene = gene, receptor_log2fc = log2FoldChange, receptor_padj = padj)

axis_fib <- lr_axis %>%
  filter(sender == "fibroblast_stromal", receiver %in% focus_receivers, ligand %in% ligand_consensus$gene, delta_healer_minus_nonhealer > 0)

axis_list <- lapply(seq_len(nrow(axis_fib)), function(i) {
  row <- axis_fib[i, ]
  rec_genes <- split_components(row$receptor_components)
  rec_de <- receptor_de_lookup %>%
    filter(receiver == row$receiver, receptor_gene %in% rec_genes)
  data.frame(
    lr_index = row$lr_index,
    interaction = row$interaction,
    ligand = row$ligand,
    receptor = row$receptor,
    receptor_components = row$receptor_components,
    receiver = row$receiver,
    curation_effort = row$curation_effort,
    mean_healer = row$mean_healer,
    mean_nonhealer = row$mean_nonhealer,
    delta_healer_minus_nonhealer = row$delta_healer_minus_nonhealer,
    ratio_healer_over_nonhealer = row$ratio_healer_over_nonhealer,
    p_wilcox = row$p_wilcox,
    padj_candidate_bh = row$padj_candidate_bh,
    receptor_component_n = length(rec_genes),
    receptor_components_with_de = collapse_unique(rec_de$receptor_gene),
    receptor_healer_high_component_n = sum(rec_de$receptor_log2fc > 0 & rec_de$receptor_padj < 0.1, na.rm = TRUE),
    receptor_mean_positive_log2fc = ifelse(nrow(rec_de) == 0, NA_real_, mean(pmax(rec_de$receptor_log2fc, 0), na.rm = TRUE))
  )
})
axis_consensus <- if (length(axis_list) == 0) data.frame() else do.call(rbind, axis_list)

axis_consensus <- axis_consensus %>%
  left_join(ligand_consensus %>% select(gene, ligand_consensus_score = consensus_score, ligand_evidence_count = evidence_count, ligand_consensus_tier = consensus_tier), by = c("ligand" = "gene")) %>%
  left_join(virtual %>% select(ligand, receiver, receiver_specific_virtual_priority = perturbation_priority, receiver_specific_knockdown_shift = predicted_knockdown_shift_toward_nonhealer), by = c("ligand", "receiver")) %>%
  mutate(
    lr_delta_component = rescale01(delta_healer_minus_nonhealer),
    lr_p_component = rescale01(pmin(-log10(pmax(p_wilcox, 1e-300)), 10)),
    ligand_component = ligand_consensus_score / 100,
    receptor_component = rescale01(receptor_mean_positive_log2fc),
    virtual_receiver_component = rescale01(receiver_specific_virtual_priority),
    curation_component = rescale01(curation_effort),
    axis_consensus_score = 100 * (
      1.4 * lr_delta_component +
        1.1 * lr_p_component +
        1.5 * ligand_component +
        0.7 * receptor_component +
        0.8 * virtual_receiver_component +
        0.5 * curation_component
    ) / 6.0
  ) %>%
  arrange(desc(axis_consensus_score), p_wilcox)

axis_out <- file.path(out_dir, "GSE165816_fibroblast_communication_axis_consensus.tsv")
write.table(axis_consensus, axis_out, sep = "\t", quote = FALSE, row.names = FALSE)

receiver_summary <- axis_consensus %>%
  group_by(receiver) %>%
  summarise(
    n_candidate_axes = n(),
    n_raw_p05_axes = sum(p_wilcox < 0.05, na.rm = TRUE),
    mean_axis_consensus_score = mean(axis_consensus_score, na.rm = TRUE),
    max_axis_consensus_score = max(axis_consensus_score, na.rm = TRUE),
    top_ligands = collapse_unique(head(unique(ligand[order(-axis_consensus_score)]), 8)),
    top_interactions = collapse_unique(head(interaction[order(-axis_consensus_score)], 8)),
    .groups = "drop"
  ) %>%
  arrange(desc(max_axis_consensus_score), desc(n_raw_p05_axes))

receiver_out <- file.path(out_dir, "GSE165816_fibroblast_communication_receiver_summary.tsv")
write.table(receiver_summary, receiver_out, sep = "\t", quote = FALSE, row.names = FALSE)

top_axes <- axis_consensus %>%
  group_by(receiver) %>%
  slice_max(axis_consensus_score, n = 6, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(axis_consensus_score)) %>%
  head(30)

if (nrow(top_axes) > 0) {
  top_axes$interaction <- factor(top_axes$interaction, levels = rev(unique(top_axes$interaction)))
  p2 <- ggplot(top_axes, aes(x = receiver, y = interaction)) +
    geom_point(aes(size = axis_consensus_score, color = ligand_consensus_tier), alpha = 0.9) +
    scale_color_manual(values = c(
      "Tier 1: convergent" = "#0B6E4F",
      "Tier 2: strong" = "#2E6F95",
      "Tier 3: supportive" = "#B75D22",
      "Exploratory" = "#777777"
    )) +
    labs(x = NULL, y = NULL, size = "Axis score", color = "Ligand tier", title = "Top fibroblast-sender communication axes") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  ggsave(file.path(fig_dir, "GSE165816_fibroblast_top_communication_axes_dotplot.png"), p2, width = 8.5, height = max(4.5, 0.18 * nrow(top_axes) + 2), dpi = 300)
}

ligand_receiver <- axis_consensus %>%
  group_by(ligand, receiver) %>%
  summarise(best_axis_score = max(axis_consensus_score, na.rm = TRUE), .groups = "drop") %>%
  left_join(ligand_consensus %>% select(gene, consensus_score), by = c("ligand" = "gene")) %>%
  arrange(desc(consensus_score), ligand)

lr_heatmap_out <- file.path(out_dir, "GSE165816_fibroblast_ligand_receiver_best_axis_scores.tsv")
write.table(ligand_receiver, lr_heatmap_out, sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(ligand_receiver) > 0) {
  ligand_levels <- ligand_consensus$gene[ligand_consensus$gene %in% ligand_receiver$ligand]
  ligand_receiver$ligand <- factor(ligand_receiver$ligand, levels = rev(ligand_levels))
  p3 <- ggplot(ligand_receiver, aes(x = receiver, y = ligand, fill = best_axis_score)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient(low = "#F2F4F3", high = "#9A3412") +
    labs(x = NULL, y = NULL, fill = "Best axis score", title = "Best fibroblast-sender LR axis by receiver") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  ggsave(file.path(fig_dir, "GSE165816_fibroblast_ligand_receiver_axis_heatmap.png"), p3, width = 7, height = max(4, 0.25 * length(ligand_levels) + 1.5), dpi = 300)
}

message("Wrote: ", ligand_out)
message("Wrote: ", evidence_out)
message("Wrote: ", axis_out)
message("Wrote: ", receiver_out)
message("Wrote figures to: ", fig_dir)
