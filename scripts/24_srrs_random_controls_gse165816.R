options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
fib_file <- file.path(project_dir, "results", "fibroblast_trajectory_gse165816", "GSE165816_fibroblast_trajectory_seurat.rds")
gene_set_file <- file.path(project_dir, "results", "srrs_framework", "SRRS_locked_gene_sets.tsv")
out_dir <- file.path(project_dir, "results", "srrs_framework")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_cells <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < 3) return(rep(NA_real_, ncol(expr)))
  mat <- as.matrix(expr[genes, , drop = FALSE])
  keep <- apply(mat, 1, sd, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 3) return(rep(NA_real_, ncol(expr)))
  colMeans(t(scale(t(mat))), na.rm = TRUE)
}

delta_by_sample <- function(cell_scores, meta) {
  df <- data.frame(
    sample_code = meta$sample_code,
    healing_status = ifelse(meta$disease == "DFU-healer", "Healer", "Non-healer"),
    score = cell_scores
  ) %>%
    filter(healing_status %in% c("Healer", "Non-healer")) %>%
    group_by(sample_code, healing_status) %>%
    summarise(sample_score = mean(score, na.rm = TRUE), .groups = "drop")
  mean(df$sample_score[df$healing_status == "Healer"], na.rm = TRUE) -
    mean(df$sample_score[df$healing_status == "Non-healer"], na.rm = TRUE)
}

message("Reading fibroblast object for random controls...")
fib <- readRDS(fib_file)
DefaultAssay(fib) <- "RNA"
expr <- LayerData(fib, assay = "RNA", layer = "data")
meta <- fib@meta.data
keep_cells <- meta$disease %in% c("DFU-healer", "DFU-nonhealer")
expr <- expr[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, ]

gene_sets <- read.delim(gene_set_file, check.names = FALSE)
sets <- split(gene_sets$gene, gene_sets$module)

detected_frac <- Matrix::rowMeans(expr > 0)
bad <- grepl("^MT-|^RPL|^RPS|^MALAT1$|^HBB|^HBA", rownames(expr))
candidate_pool <- rownames(expr)[detected_frac >= 0.05 & !bad]

target_modules <- c("stromal_ligand_panel", "d7_acute_wound_alignment", "fibroblast_repair_activation")
set.seed(20260522)
rows <- list()
random_rows <- list()

for (module in target_modules) {
  genes <- intersect(sets[[module]], rownames(expr))
  observed_delta <- delta_by_sample(score_cells(expr, genes), meta)
  n_genes <- length(genes)
  random_delta <- replicate(1000, {
    rg <- sample(setdiff(candidate_pool, genes), n_genes, replace = FALSE)
    delta_by_sample(score_cells(expr, rg), meta)
  })
  empirical_p <- (sum(random_delta >= observed_delta, na.rm = TRUE) + 1) / (sum(is.finite(random_delta)) + 1)
  rows[[module]] <- data.frame(
    module = module,
    observed_delta_healer_minus_nonhealer = observed_delta,
    random_panel_n = sum(is.finite(random_delta)),
    random_mean_delta = mean(random_delta, na.rm = TRUE),
    random_sd_delta = sd(random_delta, na.rm = TRUE),
    empirical_p_random_ge_observed = empirical_p,
    n_genes = n_genes
  )
  random_rows[[module]] <- data.frame(module = module, iteration = seq_along(random_delta), random_delta = random_delta)
}

summary_tbl <- bind_rows(rows)
random_tbl <- bind_rows(random_rows)

write.table(summary_tbl, file.path(out_dir, "GSE165816_SRRS_random_gene_panel_controls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(random_tbl, file.path(out_dir, "GSE165816_SRRS_random_gene_panel_null_distributions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote GSE165816 SRRS random controls.")
