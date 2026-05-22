options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_dir, "results", "virtual_perturbation_gse165816")
fig_dir <- file.path(project_dir, "figures", "virtual_perturbation_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

de_file <- file.path(project_dir, "results", "pseudobulk_de_gse165816", "GSE165816_all_focus_celltypes_DESeq2_healer_vs_nonhealer.tsv")
candidate_file <- file.path(project_dir, "results", "integrated_candidates_gse165816", "GSE165816_integrated_fibroblast_sender_ligand_shortlist.tsv")
ltm_file <- file.path(project_dir, "data", "metadata", "nichenet_ligand_target_matrix_zenodo3260758.rds")

de <- read.delim(de_file, check.names = FALSE)
candidates <- read.delim(candidate_file, check.names = FALSE)
ligand_target_matrix <- readRDS(ltm_file)

program_sets <- list(
  fibroblast_sender_ligand_shortlist = c("IL11", "CCL20", "INHBA", "SERPINE1", "IL6", "THBS1", "TNC", "WNT5A", "ADAM12", "PTGS2"),
  fibroblast_repair_activation = c("THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20", "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A", "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2"),
  ecm_remodeling_migration = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1", "COL6A2", "COL10A1", "THBS1", "TNC", "MMP1", "MMP3", "MMP10", "SERPINE1", "ITGA5", "ITGAV", "PDPN", "TNFAIP6", "CTHRC1"),
  repair_inflammatory_signaling = c("IL6", "IL11", "INHBA", "AIM2", "CCL20", "CXCL5", "CXCL8", "TNFAIP3", "TNFRSF12A", "PTGS2", "FPR2", "TNIP3", "SLC39A8"),
  resolution_metabolic_myeloid = c("FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CES1", "CXCL5", "TNC", "COL4A1", "COL4A2")
)

priority_ligands <- c("IL11", "CCL20", "INHBA", "SERPINE1", "IL6", "THBS1", "TNC", "WNT5A", "ADAM12", "PTGS2", "TGFB1", "FGF2", "VEGFA", "HGF")
priority_ligands <- intersect(priority_ligands, colnames(ligand_target_matrix))

de_unique <- de %>%
  filter(is.finite(log2FoldChange)) %>%
  mutate(
    padj_filled = ifelse(is.na(padj), 1, padj),
    pvalue_filled = ifelse(is.na(pvalue), 1, pvalue)
  ) %>%
  arrange(cell_type, gene, padj_filled, pvalue_filled, desc(abs(log2FoldChange))) %>%
  group_by(cell_type, gene) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-padj_filled, -pvalue_filled)

ligand_target_tbl <- bind_rows(lapply(priority_ligands, function(lig) {
  data.frame(
    ligand = lig,
    target = rownames(ligand_target_matrix),
    weight = as.numeric(ligand_target_matrix[, lig])
  )
})) %>%
  filter(is.finite(weight), weight > 0) %>%
  group_by(ligand, target) %>%
  summarise(weight = max(weight, na.rm = TRUE), .groups = "drop")

receivers <- sort(unique(de_unique$cell_type))

links <- bind_rows(lapply(priority_ligands, function(lig) {
  bind_rows(lapply(receivers, function(receiver) {
    dz <- de_unique %>%
      filter(cell_type == receiver, is.finite(log2FoldChange)) %>%
      arrange(desc(log2FoldChange))
    z <- ligand_target_tbl %>%
      filter(ligand == lig, target %in% dz$gene) %>%
      arrange(desc(weight)) %>%
      slice_head(n = 250) %>%
      select(target, weight)
    if (nrow(z) == 0) {
      return(data.frame())
    }
    z$ligand <- lig
    z$sender <- "fibroblast_stromal"
    z$receiver <- receiver
    z
  }))
}))

score_pair <- function(ligand_id, receiver_id) {
  z <- links %>% filter(.data$ligand == ligand_id, .data$receiver == receiver_id)
  if (nrow(z) == 0) {
    return(data.frame())
  }
  dz <- de_unique %>% filter(.data$cell_type == receiver_id) %>% select(gene, log2FoldChange, pvalue, padj, direction)
  zz <- z %>% inner_join(dz, by = c("target" = "gene"))
  if (nrow(zz) == 0) {
    return(data.frame())
  }
  zz$padj_filled <- ifelse(is.na(zz$padj), 1, zz$padj)
  weighted_lfc <- sum(zz$weight * zz$log2FoldChange, na.rm = TRUE) / sum(zz$weight, na.rm = TRUE)
  weighted_abs_lfc <- sum(zz$weight * abs(zz$log2FoldChange), na.rm = TRUE) / sum(zz$weight, na.rm = TRUE)
  healer_high_weight_fraction <- sum(zz$weight[zz$log2FoldChange > 0], na.rm = TRUE) / sum(zz$weight, na.rm = TRUE)
  healer_high_sig_weight_fraction <- sum(zz$weight[zz$log2FoldChange > 0 & zz$padj_filled < 0.1], na.rm = TRUE) / sum(zz$weight, na.rm = TRUE)
  weighted_pseudo_impact <- weighted_lfc * log2(1 + nrow(zz))

  program_rows <- lapply(names(program_sets), function(pn) {
    genes <- program_sets[[pn]]
    ov <- zz %>% filter(target %in% genes)
    data.frame(
      ligand = ligand_id,
      receiver = receiver_id,
      program = pn,
      target_overlap_n = nrow(ov),
      target_overlap_genes = paste(unique(ov$target), collapse = ";"),
      target_overlap_weight = sum(ov$weight, na.rm = TRUE),
      target_overlap_weighted_lfc = ifelse(nrow(ov) > 0, sum(ov$weight * ov$log2FoldChange, na.rm = TRUE) / sum(ov$weight, na.rm = TRUE), NA_real_)
    )
  })
  program_df <- bind_rows(program_rows)

  data.frame(
    ligand = ligand_id,
    receiver = receiver_id,
    target_n = nrow(zz),
    target_genes = paste(unique(zz$target), collapse = ";"),
    top_targets = paste(head(zz$target[order(-zz$weight)], 12), collapse = ";"),
    weighted_target_log2fc_healer_vs_nonhealer = weighted_lfc,
    weighted_abs_target_log2fc = weighted_abs_lfc,
    healer_high_target_weight_fraction = healer_high_weight_fraction,
    healer_high_sig_target_weight_fraction = healer_high_sig_weight_fraction,
    predicted_knockdown_shift_toward_nonhealer = weighted_pseudo_impact,
    program_overlap_summary = paste(program_df$program[program_df$target_overlap_n > 0], program_df$target_overlap_genes[program_df$target_overlap_n > 0], sep = ":", collapse = " | ")
  )
}

summary <- bind_rows(lapply(priority_ligands, function(ligand_id) {
  bind_rows(lapply(sort(unique(links$receiver)), function(receiver_id) score_pair(ligand_id, receiver_id)))
}))

summary <- summary %>%
  left_join(
    candidates %>%
      select(gene, priority_score, log2FoldChange, padj, nichenet_max_pearson, nichenet_best_receiver, bulk_delta_healer_minus_nonhealer),
    by = c("ligand" = "gene")
  ) %>%
  mutate(
    perturbation_priority = predicted_knockdown_shift_toward_nonhealer *
      ifelse(is.na(priority_score), 1, priority_score) *
      ifelse(healer_high_sig_target_weight_fraction > 0, 1.25, 1)
  ) %>%
  arrange(desc(perturbation_priority))

write.table(
  summary,
  file.path(out_dir, "GSE165816_virtual_ligand_knockdown_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

program_details <- bind_rows(lapply(priority_ligands, function(ligand_id) {
  bind_rows(lapply(sort(unique(links$receiver)), function(receiver_id) {
    z <- links %>% filter(.data$ligand == ligand_id, .data$receiver == receiver_id)
    if (nrow(z) == 0) return(data.frame())
    dz <- de_unique %>% filter(.data$cell_type == receiver_id) %>% select(gene, log2FoldChange, pvalue, padj, direction)
    zz <- z %>% inner_join(dz, by = c("target" = "gene"))
    if (nrow(zz) == 0) return(data.frame())
    bind_rows(lapply(names(program_sets), function(pn) {
      ov <- zz %>% filter(target %in% program_sets[[pn]])
      data.frame(
        ligand = ligand_id,
        receiver = receiver_id,
        program = pn,
        overlap_n = nrow(ov),
        overlap_genes = paste(unique(ov$target), collapse = ";"),
        overlap_weight = sum(ov$weight, na.rm = TRUE),
        weighted_log2fc = ifelse(nrow(ov) > 0, sum(ov$weight * ov$log2FoldChange, na.rm = TRUE) / sum(ov$weight, na.rm = TRUE), NA_real_)
      )
    }))
  }))
}))

write.table(
  program_details,
  file.path(out_dir, "GSE165816_virtual_ligand_knockdown_program_overlaps.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top <- summary %>% filter(target_n >= 3) %>% slice_head(n = 30)
if (nrow(top) > 0) {
  top$label <- paste(top$ligand, "->", top$receiver)
  top$label <- factor(top$label, levels = rev(top$label))
  p <- ggplot(top, aes(x = label, y = perturbation_priority, fill = receiver)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Predicted knockdown shift toward nonhealer", fill = "Receiver") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(fig_dir, "GSE165816_virtual_ligand_knockdown_priority.png"),
    p,
    width = 8.5,
    height = 7,
    dpi = 220
  )
}

heat <- summary %>%
  filter(target_n >= 3) %>%
  mutate(label = paste(ligand, receiver, sep = "->")) %>%
  select(ligand, receiver, predicted_knockdown_shift_toward_nonhealer)
if (nrow(heat) > 0) {
  p2 <- ggplot(heat, aes(x = receiver, y = ligand, fill = predicted_knockdown_shift_toward_nonhealer)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient2(low = "#D55E00", mid = "#F5F5F5", high = "#0072B2", midpoint = 0) +
    labs(x = "Receiver", y = "Fibroblast ligand", fill = "Predicted\nKD shift") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(
    file.path(fig_dir, "GSE165816_virtual_ligand_knockdown_heatmap.png"),
    p2,
    width = 7,
    height = 5.5,
    dpi = 220
  )
}

message("Wrote virtual perturbation outputs.")
