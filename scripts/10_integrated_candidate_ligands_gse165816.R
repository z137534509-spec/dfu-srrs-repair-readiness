options(stringsAsFactors = FALSE)

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_dir, "results", "integrated_candidates_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fib_de_file <- file.path(project_dir, "results", "pseudobulk_de_gse165816", "GSE165816_fibroblast_stromal_DESeq2_healer_vs_nonhealer.tsv")
nichenet_file <- file.path(project_dir, "results", "nichenet_gse165816", "GSE165816_nichenet_ligand_activities.tsv")
lr_file <- file.path(project_dir, "results", "omnipath_lr_gse165816", "GSE165816_omnipath_lr_healer_vs_nonhealer_tests.tsv")
bulk_file <- file.path(project_dir, "results", "bulk_validation_gse134431", "GSE134431_candidate_gene_healer_vs_nonhealer_tests.tsv")

fib_de <- read.delim(fib_de_file, check.names = FALSE)
nichenet <- read.delim(nichenet_file, check.names = FALSE)
lr <- read.delim(lr_file, check.names = FALSE)
bulk <- read.delim(bulk_file, check.names = FALSE)

fib_niche <- nichenet[nichenet$sender == "fibroblast_stromal", ]
fib_niche <- fib_niche[fib_niche$receiver %in% c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte"), ]

summarize_ligand <- function(lig) {
  z <- fib_niche[fib_niche$test_ligand == lig, ]
  z <- z[order(-z$pearson), ]
  data.frame(
    gene = lig,
    nichenet_max_pearson = max(z$pearson, na.rm = TRUE),
    nichenet_mean_pearson = mean(z$pearson, na.rm = TRUE),
    nichenet_best_receiver = z$receiver[1],
    nichenet_best_auroc = z$auroc[1],
    nichenet_receivers_positive = paste(z$receiver[z$pearson > 0], collapse = ";")
  )
}
niche_summary <- do.call(rbind, lapply(sort(unique(fib_niche$test_ligand)), summarize_ligand))

lr_fib <- lr[lr$sender == "fibroblast_stromal" & lr$delta_healer_minus_nonhealer > 0, ]
lr_summary <- do.call(rbind, lapply(sort(unique(lr_fib$ligand)), function(lig) {
  z <- lr_fib[lr_fib$ligand == lig, ]
  z <- z[order(z$p_wilcox, -z$ratio_healer_over_nonhealer), ]
  data.frame(
    gene = lig,
    lr_best_interaction = z$interaction[1],
    lr_best_receiver = z$receiver[1],
    lr_best_raw_p = z$p_wilcox[1],
    lr_best_ratio = z$ratio_healer_over_nonhealer[1],
    lr_healer_higher_pairs_raw_p05 = sum(z$p_wilcox < 0.05, na.rm = TRUE)
  )
}))

merged <- merge(niche_summary, fib_de[, c("gene", "log2FoldChange", "pvalue", "padj", "direction")], by = "gene", all.x = TRUE)
merged <- merge(merged, lr_summary, by = "gene", all.x = TRUE)
merged <- merge(
  merged,
  bulk[, c("gene", "delta_healer_minus_nonhealer", "p_wilcox", "padj_bh")],
  by = "gene",
  all.x = TRUE,
  suffixes = c("", "_bulk")
)
names(merged)[names(merged) == "delta_healer_minus_nonhealer"] <- "bulk_delta_healer_minus_nonhealer"
names(merged)[names(merged) == "p_wilcox"] <- "bulk_p_wilcox"
names(merged)[names(merged) == "padj_bh"] <- "bulk_padj_bh"

merged$fibroblast_healer_high <- !is.na(merged$log2FoldChange) & merged$log2FoldChange > 0 & !is.na(merged$padj) & merged$padj < 0.1
merged$nichenet_support <- !is.na(merged$nichenet_max_pearson) & merged$nichenet_max_pearson > 0.04
merged$lr_support <- !is.na(merged$lr_best_raw_p) & merged$lr_best_raw_p < 0.05
merged$bulk_direction_support <- !is.na(merged$bulk_delta_healer_minus_nonhealer) & merged$bulk_delta_healer_minus_nonhealer > 0

merged$priority_score <- 0
merged$priority_score <- merged$priority_score + ifelse(merged$fibroblast_healer_high, 3, 0)
merged$priority_score <- merged$priority_score + ifelse(merged$nichenet_support, 2, 0)
merged$priority_score <- merged$priority_score + ifelse(merged$lr_support, 2, 0)
merged$priority_score <- merged$priority_score + ifelse(merged$bulk_direction_support, 1, 0)
merged$priority_score <- merged$priority_score + ifelse(!is.na(merged$padj) & merged$padj < 0.01, 1, 0)

merged <- merged[order(-merged$priority_score, -merged$nichenet_max_pearson, merged$padj), ]

out_file <- file.path(out_dir, "GSE165816_integrated_fibroblast_sender_ligand_priority.tsv")
write.table(merged, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

shortlist <- merged[merged$priority_score >= 6, ]
shortlist <- shortlist[
  ,
  c(
    "gene", "priority_score", "log2FoldChange", "padj", "nichenet_max_pearson",
    "nichenet_best_receiver", "nichenet_best_auroc", "lr_best_interaction",
    "lr_best_receiver", "lr_best_raw_p", "lr_best_ratio",
    "bulk_delta_healer_minus_nonhealer", "bulk_p_wilcox",
    "fibroblast_healer_high", "nichenet_support", "lr_support", "bulk_direction_support"
  )
]
shortlist_file <- file.path(out_dir, "GSE165816_integrated_fibroblast_sender_ligand_shortlist.tsv")
write.table(shortlist, shortlist_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", out_file)
message("Wrote: ", shortlist_file)
