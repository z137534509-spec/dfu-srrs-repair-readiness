options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
infile <- file.path(
  project_dir,
  "data", "raw", "GSE134431",
  "GSE134431_181016_diabetic_rna-seq_results.gene.rpkm.xlsx"
)
out_dir <- file.path(project_dir, "results", "bulk_validation_gse134431")
fig_dir <- file.path(project_dir, "figures", "bulk_validation_gse134431")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

raw <- readxl::read_excel(infile, sheet = "181016_diabetic_rna-seq_results", col_names = FALSE)
sample_ids <- as.character(unlist(raw[3, ]))
sample_cols <- which(grepl("^s[0-9]+$", sample_ids))
genes <- toupper(as.character(raw[[1]][-(1:3)]))
expr <- raw[-(1:3), sample_cols]
expr <- as.data.frame(lapply(expr, function(v) suppressWarnings(as.numeric(v))))
names(expr) <- sample_ids[sample_cols]
expr$gene <- genes
expr <- expr[!is.na(expr$gene) & expr$gene != "", ]
expr <- aggregate(. ~ gene, expr, mean, na.rm = TRUE)
rownames(expr) <- expr$gene
expr$gene <- NULL
expr <- as.matrix(expr)
expr_log <- log2(expr + 1)

samples <- as.data.frame(readxl::read_excel(infile, sheet = "Samples"))
names(samples) <- make.names(names(samples))
samples$sample_id <- samples$Sample.ID
samples$validation_group <- ifelse(samples$Group == "DFS", "DFS", NA_character_)
samples$validation_group[samples$sample_id %in% paste0("s", 14:20)] <- "DFU_healer"
samples$validation_group[samples$sample_id %in% c(paste0("s", 21:25), "s27")] <- "DFU_nonhealer"
samples <- samples[samples$sample_id %in% colnames(expr_log), ]

signature_sets <- list(
  fibroblast_repair_activation = c(
    "THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20",
    "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A",
    "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2"
  ),
  ecm_remodeling_migration = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
    "COL6A1", "COL6A2", "COL10A1", "THBS1", "TNC", "MMP1", "MMP3",
    "MMP10", "SERPINE1", "ITGA5", "ITGAV", "PDPN", "TNFAIP6", "CTHRC1"
  ),
  repair_inflammatory_signaling = c(
    "IL6", "IL11", "INHBA", "AIM2", "CCL20", "CXCL5", "CXCL8",
    "TNFAIP3", "TNFRSF12A", "PTGS2", "FPR2", "TNIP3", "SLC39A8"
  ),
  resolution_metabolic_myeloid = c(
    "FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8",
    "CES1", "CXCL5", "TNC", "COL4A1", "COL4A2"
  )
)

score_signature <- function(signature_name, geneset) {
  present <- intersect(toupper(geneset), rownames(expr_log))
  m <- expr_log[present, , drop = FALSE]
  z <- t(apply(m, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      rep(0, length(v))
    } else {
      (v - mean(v, na.rm = TRUE)) / s
    }
  }))
  data.frame(
    sample_id = colnames(expr_log),
    signature = signature_name,
    score = colMeans(z, na.rm = TRUE),
    genes_present = length(present),
    genes_present_list = paste(present, collapse = ";")
  )
}

scores <- do.call(rbind, Map(score_signature, names(signature_sets), signature_sets))
scores <- merge(scores, samples, by = "sample_id", all.x = TRUE)

scores_out <- file.path(out_dir, "GSE134431_signature_scores.tsv")
write.table(scores, scores_out, sep = "\t", quote = FALSE, row.names = FALSE)

dfu_scores <- scores[scores$validation_group %in% c("DFU_healer", "DFU_nonhealer"), ]
test_one <- function(sig) {
  z <- dfu_scores[dfu_scores$signature == sig, ]
  healer <- z$score[z$validation_group == "DFU_healer"]
  nonhealer <- z$score[z$validation_group == "DFU_nonhealer"]
  data.frame(
    signature = sig,
    genes_present = unique(z$genes_present)[1],
    mean_healer = mean(healer),
    mean_nonhealer = mean(nonhealer),
    delta_healer_minus_nonhealer = mean(healer) - mean(nonhealer),
    p_wilcox = tryCatch(wilcox.test(healer, nonhealer, exact = FALSE)$p.value, error = function(e) NA_real_),
    p_ttest = tryCatch(t.test(healer, nonhealer)$p.value, error = function(e) NA_real_)
  )
}
sig_tests <- do.call(rbind, lapply(names(signature_sets), test_one))
sig_tests$padj_wilcox_bh <- p.adjust(sig_tests$p_wilcox, method = "BH")
sig_tests <- sig_tests[order(sig_tests$p_wilcox), ]

sig_tests_out <- file.path(out_dir, "GSE134431_signature_healer_vs_nonhealer_tests.tsv")
write.table(sig_tests, sig_tests_out, sep = "\t", quote = FALSE, row.names = FALSE)

candidate_genes <- unique(toupper(unlist(signature_sets)))
candidate_genes <- intersect(candidate_genes, rownames(expr_log))
gene_tests <- do.call(rbind, lapply(candidate_genes, function(g) {
  vals <- data.frame(sample_id = colnames(expr_log), gene = g, expr_log2_rpkm1 = as.numeric(expr_log[g, ]))
  vals <- merge(vals, samples, by = "sample_id")
  vals <- vals[vals$validation_group %in% c("DFU_healer", "DFU_nonhealer"), ]
  healer <- vals$expr_log2_rpkm1[vals$validation_group == "DFU_healer"]
  nonhealer <- vals$expr_log2_rpkm1[vals$validation_group == "DFU_nonhealer"]
  data.frame(
    gene = g,
    mean_healer = mean(healer),
    mean_nonhealer = mean(nonhealer),
    delta_healer_minus_nonhealer = mean(healer) - mean(nonhealer),
    p_wilcox = tryCatch(wilcox.test(healer, nonhealer, exact = FALSE)$p.value, error = function(e) NA_real_)
  )
}))
gene_tests$padj_bh <- p.adjust(gene_tests$p_wilcox, method = "BH")
gene_tests <- gene_tests[order(gene_tests$p_wilcox), ]

gene_tests_out <- file.path(out_dir, "GSE134431_candidate_gene_healer_vs_nonhealer_tests.tsv")
write.table(gene_tests, gene_tests_out, sep = "\t", quote = FALSE, row.names = FALSE)

plot_scores <- scores[scores$validation_group %in% c("DFS", "DFU_healer", "DFU_nonhealer"), ]
plot_scores$validation_group <- factor(plot_scores$validation_group, levels = c("DFS", "DFU_healer", "DFU_nonhealer"))
plot_scores$signature <- factor(plot_scores$signature, levels = names(signature_sets))

p <- ggplot(plot_scores, aes(x = validation_group, y = score, fill = validation_group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.75) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.8, alpha = 0.8) +
  facet_wrap(~ signature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(DFS = "#7F7F7F", DFU_healer = "#0072B2", DFU_nonhealer = "#D55E00")) +
  labs(x = NULL, y = "Mean gene-wise z-score", fill = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(
  file.path(fig_dir, "GSE134431_signature_scores_boxplot.png"),
  p,
  width = 8,
  height = 5.2,
  dpi = 220
)

top_genes <- head(gene_tests, 30)
top_genes$gene <- factor(top_genes$gene, levels = rev(top_genes$gene))
p2 <- ggplot(top_genes, aes(x = gene, y = delta_healer_minus_nonhealer, fill = delta_healer_minus_nonhealer > 0)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00")) +
  labs(x = NULL, y = "Delta log2(RPKM + 1), healer - nonhealer", fill = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")
ggsave(
  file.path(fig_dir, "GSE134431_top_candidate_gene_deltas.png"),
  p2,
  width = 7,
  height = 6,
  dpi = 220
)

message("Wrote: ", scores_out)
message("Wrote: ", sig_tests_out)
message("Wrote: ", gene_tests_out)
