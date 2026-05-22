options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(jsonlite)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "raw", "GSE166120")
dcc_dir <- file.path(raw_dir, "dcc")
pkc_file <- file.path(raw_dir, "GSE166120_Alpha_CTA_v2.0.pkc.txt.gz")
series_file <- file.path(project_dir, "data", "metadata", "GSE166120_series_matrix.txt.gz")

out_dir <- file.path(project_dir, "results", "geomx_validation_gse166120")
fig_dir <- file.path(project_dir, "figures", "geomx_validation_gse166120")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

strip_quotes <- function(x) {
  gsub('^"|"$', "", x)
}

parse_series_line <- function(lines, prefix) {
  line <- grep(paste0("^", prefix), lines, value = TRUE)
  if (length(line) == 0) return(character())
  strip_quotes(strsplit(line[1], "\t", fixed = TRUE)[[1]][-1])
}

series_lines <- readLines(gzfile(series_file), warn = FALSE)
sample_ids <- parse_series_line(series_lines, "!Sample_geo_accession")
sample_titles <- parse_series_line(series_lines, "!Sample_title")
meta <- data.frame(sample_id = sample_ids, roi_id = sample_titles)

char_lines <- grep("^!Sample_characteristics_ch1", series_lines, value = TRUE)
for (line in char_lines) {
  vals <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
  keys <- sub(":.*$", "", vals)
  values <- sub("^[^:]+:\\s*", "", vals)
  key <- tolower(gsub("[^A-Za-z0-9]+", "_", keys[1]))
  meta[[key]] <- values
}

meta <- meta %>%
  mutate(
    healing_status = factor(healing_status, levels = c("Non-Healer", "Healer")),
    roi_location_clean = tolower(gsub("[ -]+", "_", roi_location)),
    wound_roi = !grepl("non_ulcer|epidermis", roi_location_clean)
  )

pkc <- jsonlite::fromJSON(pkc_file, simplifyDataFrame = FALSE)
probe_map <- bind_rows(lapply(pkc[["Targets"]], function(target) {
  probes <- target[["Probes"]]
  if (length(probes) == 0) return(data.frame())
  data.frame(
    gene = target[["DisplayName"]],
    code_class = target[["CodeClass"]],
    rts_id = vapply(probes, function(probe) probe[["RTS_ID"]], character(1))
  )
})) %>%
  filter(grepl("^Endogenous", code_class), !is.na(rts_id), rts_id != "") %>%
  distinct(rts_id, gene, .keep_all = TRUE)

parse_dcc <- function(file) {
  lines <- readLines(gzfile(file), warn = FALSE)
  start <- grep("^<Code_Summary>", lines)
  end <- grep("^</Code_Summary>", lines)
  if (length(start) == 0 || length(end) == 0 || end <= start) {
    stop("No Code_Summary section in ", file)
  }
  tab <- read.csv(
    text = paste(lines[(start + 1):(end - 1)], collapse = "\n"),
    header = FALSE,
    col.names = c("rts_id", "count")
  )
  sample_id <- sub("_.*$", "", basename(file))
  tab$sample_id <- sample_id
  tab
}

dcc_files <- list.files(dcc_dir, pattern = "\\.dcc\\.gz$", full.names = TRUE)
probe_counts <- bind_rows(lapply(dcc_files, parse_dcc)) %>%
  inner_join(probe_map, by = "rts_id") %>%
  group_by(gene, sample_id) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

count_mat <- probe_counts %>%
  tidyr::pivot_wider(names_from = sample_id, values_from = count, values_fill = 0) %>%
  as.data.frame()
rownames(count_mat) <- count_mat$gene
count_mat$gene <- NULL
count_mat <- as.matrix(count_mat)
count_mat <- count_mat[, sample_ids, drop = FALSE]

uq <- apply(count_mat, 2, function(x) as.numeric(quantile(x[x > 0], 0.75, na.rm = TRUE)))
uq[!is.finite(uq) | uq <= 0] <- median(uq[uq > 0], na.rm = TRUE)
norm_mat <- t(t(count_mat) / uq * median(uq, na.rm = TRUE))
log_mat <- log2(norm_mat + 1)

program_sets <- list(
  fibroblast_sender_ligand_shortlist = c("IL11", "CCL20", "INHBA", "SERPINE1", "IL6", "THBS1", "TNC", "WNT5A", "ADAM12", "PTGS2"),
  fibroblast_repair_activation = c("THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20", "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A", "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2"),
  ecm_remodeling_migration = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1", "COL6A2", "COL10A1", "THBS1", "TNC", "MMP1", "MMP3", "MMP10", "SERPINE1", "ITGA5", "ITGAV", "PDPN", "TNFAIP6", "CTHRC1"),
  repair_inflammatory_signaling = c("IL6", "IL11", "INHBA", "AIM2", "CCL20", "CXCL5", "CXCL8", "TNFAIP3", "TNFRSF12A", "PTGS2", "FPR2", "TNIP3", "SLC39A8"),
  resolution_metabolic_myeloid = c("FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CES1", "CXCL5", "TNC", "COL4A1", "COL4A2"),
  published_healing_enriched_fibroblast = c("MMP1", "MMP3", "MMP11", "HIF1A", "CHI3L1", "TNFAIP6"),
  gp130_stromal = c("IL6", "IL11", "IL6ST", "STAT3", "SOCS3", "JUNB"),
  tgf_like_repair = c("TGFB1", "INHBA", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF")
)

score_signature <- function(genes) {
  genes <- intersect(genes, rownames(log_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(log_mat)))
  }
  mat <- log_mat[genes, , drop = FALSE]
  keep <- apply(mat, 1, sd, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0) {
    return(rep(NA_real_, ncol(log_mat)))
  }
  z <- t(scale(t(mat)))
  colMeans(z, na.rm = TRUE)
}

signature_coverage <- bind_rows(lapply(names(program_sets), function(sig) {
  data.frame(
    signature = sig,
    requested_n = length(program_sets[[sig]]),
    available_n = length(intersect(program_sets[[sig]], rownames(log_mat))),
    available_genes = paste(intersect(program_sets[[sig]], rownames(log_mat)), collapse = ";")
  )
}))

score_mat <- sapply(program_sets, score_signature)
score_df <- as.data.frame(score_mat) %>%
  mutate(sample_id = rownames(.)) %>%
  tidyr::pivot_longer(-sample_id, names_to = "signature", values_to = "score") %>%
  left_join(meta, by = "sample_id")

roi_stats <- score_df %>%
  filter(wound_roi, !is.na(score)) %>%
  group_by(signature) %>%
  summarise(
    healer_roi_n = sum(healing_status == "Healer"),
    nonhealer_roi_n = sum(healing_status == "Non-Healer"),
    healer_mean = mean(score[healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(score[healing_status == "Non-Healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = healer_mean - nonhealer_mean,
    p_wilcox_roi = suppressWarnings(wilcox.test(score ~ healing_status)$p.value),
    .groups = "drop"
  ) %>%
  mutate(BH_roi = p.adjust(p_wilcox_roi, method = "BH"))

subject_scores <- score_df %>%
  filter(wound_roi, !is.na(score)) %>%
  group_by(subject_id, healing_status, signature) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), roi_n = n(), .groups = "drop")

subject_stats <- subject_scores %>%
  group_by(signature) %>%
  summarise(
    healer_subject_n = sum(healing_status == "Healer"),
    nonhealer_subject_n = sum(healing_status == "Non-Healer"),
    healer_mean = mean(mean_score[healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(mean_score[healing_status == "Non-Healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = healer_mean - nonhealer_mean,
    p_wilcox_subject = suppressWarnings(wilcox.test(mean_score ~ healing_status)$p.value),
    .groups = "drop"
  ) %>%
  mutate(BH_subject = p.adjust(p_wilcox_subject, method = "BH"))

candidate_genes <- c("IL6", "IL11", "CCL20", "INHBA", "TGFB1", "TNC", "PTGS2", "THBS1", "SERPINE1", "ADAM12", "WNT5A", "VEGFA", "HGF", "FGF2")
candidate_presence <- data.frame(
  gene = candidate_genes,
  available = candidate_genes %in% rownames(log_mat)
)

candidate_expr <- log_mat[intersect(candidate_genes, rownames(log_mat)), , drop = FALSE] %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "log_expr") %>%
  left_join(meta, by = "sample_id")

candidate_subject <- candidate_expr %>%
  filter(wound_roi) %>%
  group_by(subject_id, healing_status, gene) %>%
  summarise(mean_log_expr = mean(log_expr, na.rm = TRUE), roi_n = n(), .groups = "drop")

candidate_subject_stats <- candidate_subject %>%
  group_by(gene) %>%
  summarise(
    healer_subject_n = sum(healing_status == "Healer"),
    nonhealer_subject_n = sum(healing_status == "Non-Healer"),
    healer_mean = mean(mean_log_expr[healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(mean_log_expr[healing_status == "Non-Healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = healer_mean - nonhealer_mean,
    p_wilcox_subject = suppressWarnings(wilcox.test(mean_log_expr ~ healing_status)$p.value),
    .groups = "drop"
  ) %>%
  mutate(BH_subject = p.adjust(p_wilcox_subject, method = "BH"))

write.table(meta, file.path(out_dir, "GSE166120_geomx_roi_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(signature_coverage, file.path(out_dir, "GSE166120_geomx_signature_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(score_df, file.path(out_dir, "GSE166120_geomx_signature_scores_by_roi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subject_scores, file.path(out_dir, "GSE166120_geomx_signature_scores_by_subject.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(roi_stats, file.path(out_dir, "GSE166120_geomx_signature_healer_vs_nonhealer_roi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subject_stats, file.path(out_dir, "GSE166120_geomx_signature_healer_vs_nonhealer_subject.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(candidate_presence, file.path(out_dir, "GSE166120_geomx_candidate_gene_presence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(candidate_expr, file.path(out_dir, "GSE166120_geomx_candidate_gene_expression_by_roi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(candidate_subject_stats, file.path(out_dir, "GSE166120_geomx_candidate_gene_healer_vs_nonhealer_subject.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_sigs <- c(
  "fibroblast_sender_ligand_shortlist",
  "fibroblast_repair_activation",
  "ecm_remodeling_migration",
  "published_healing_enriched_fibroblast",
  "repair_inflammatory_signaling",
  "resolution_metabolic_myeloid",
  "gp130_stromal",
  "tgf_like_repair"
)

p1 <- score_df %>%
  filter(wound_roi, signature %in% plot_sigs) %>%
  ggplot(aes(x = healing_status, y = score, color = healing_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.2) +
  geom_jitter(aes(shape = subject_id), width = 0.12, height = 0, size = 2.2, alpha = 0.9) +
  facet_wrap(~signature, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Non-Healer" = "#D55E00", "Healer" = "#0072B2")) +
  labs(x = NULL, y = "ROI signature score", color = NULL, shape = "Subject") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "GSE166120_geomx_signature_scores_ulcer_roi.png"), p1, width = 9, height = 6.5, dpi = 220)

p2 <- candidate_subject %>%
  filter(gene %in% intersect(candidate_genes, rownames(log_mat))) %>%
  ggplot(aes(x = healing_status, y = mean_log_expr, color = healing_status)) +
  geom_point(aes(shape = subject_id), size = 2.6, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  scale_color_manual(values = c("Non-Healer" = "#D55E00", "Healer" = "#0072B2")) +
  labs(x = NULL, y = "Subject mean log2 normalized expression", color = NULL, shape = "Subject") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "GSE166120_geomx_candidate_ligand_expression_subject.png"), p2, width = 9, height = 6.5, dpi = 220)

message("Wrote GSE166120 GeoMx validation outputs.")
