options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
seurat_file <- file.path(
  project_dir,
  "results", "seurat_gse165816_foot_skin",
  "GSE165816_foot_skin_seurat_firstpass.rds"
)
lr_file <- file.path(
  project_dir,
  "data", "metadata",
  "omnipath_ligrecextra_interactions_2026-05-21.tsv"
)
out_dir <- file.path(project_dir, "results", "omnipath_lr_gse165816")
fig_dir <- file.path(project_dir, "figures", "omnipath_lr_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

min_curation <- 2
focus_celltypes <- c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte", "pericyte_smc")
focus_diseases <- c("DFU-healer", "DFU-nonhealer")

split_complex <- function(x) {
  x <- gsub("^COMPLEX:", "", x)
  parts <- unlist(strsplit(x, "_", fixed = TRUE), use.names = FALSE)
  unique(parts[nzchar(parts)])
}

message("Reading OmniPath LR table...")
lr <- read.delim(lr_file, check.names = FALSE)
lr <- lr[
  !is.na(lr$source_genesymbol) &
    !is.na(lr$target_genesymbol) &
    lr$source_genesymbol != "" &
    lr$target_genesymbol != "" &
    lr$curation_effort >= min_curation,
]
lr$lr_index <- seq_len(nrow(lr))
lr$source_components <- vapply(lr$source_genesymbol, function(x) paste(split_complex(x), collapse = ";"), character(1))
lr$target_components <- vapply(lr$target_genesymbol, function(x) paste(split_complex(x), collapse = ";"), character(1))
lr$interaction <- paste(lr$source_genesymbol, lr$target_genesymbol, sep = "->")

all_components <- unique(unlist(strsplit(paste(c(lr$source_components, lr$target_components), collapse = ";"), ";", fixed = TRUE)))
all_components <- all_components[nzchar(all_components)]

message("Loading Seurat object...")
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
mat <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
)

components_present <- intersect(all_components, rownames(mat))
message("LR components present in Seurat object: ", length(components_present), "/", length(all_components))

source_ok <- vapply(strsplit(lr$source_components, ";", fixed = TRUE), function(g) all(g %in% components_present), logical(1))
target_ok <- vapply(strsplit(lr$target_components, ";", fixed = TRUE), function(g) all(g %in% components_present), logical(1))
lr$components_present <- source_ok & target_ok
lr <- lr[lr$components_present, ]
lr$lr_index <- seq_len(nrow(lr))
message("LR interactions retained after component filter: ", nrow(lr))

meta <- obj@meta.data
meta$cell_id <- rownames(meta)
meta <- meta[meta$disease %in% focus_diseases & meta$broad_cell_type %in% focus_celltypes, ]
meta$sample_id <- if ("sample_code" %in% names(meta)) meta$sample_code else meta$geo_accession

groups <- unique(meta[, c("sample_id", "geo_accession", "disease", "broad_cell_type")])
groups <- groups[order(groups$disease, groups$sample_id, groups$broad_cell_type), ]
groups$group_id <- seq_len(nrow(groups))
groups$group_key <- paste(groups$sample_id, groups$broad_cell_type, sep = "\r")

message("Summarizing LR component expression by sample and cell type...")
expr_summary_list <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  g <- groups[i, ]
  cells <- meta$cell_id[
    meta$sample_id == g$sample_id &
      meta$disease == g$disease &
      meta$broad_cell_type == g$broad_cell_type
  ]
  m <- mat[components_present, cells, drop = FALSE]
  expr_summary_list[[i]] <- data.frame(
    group_id = g$group_id,
    sample_id = g$sample_id,
    geo_accession = g$geo_accession,
    disease = g$disease,
    broad_cell_type = g$broad_cell_type,
    gene = components_present,
    avg_expr = as.numeric(Matrix::rowMeans(m)),
    pct_expr = as.numeric(Matrix::rowMeans(m > 0)),
    cells = length(cells)
  )
}
expr_summary <- do.call(rbind, expr_summary_list)

expr_out <- file.path(out_dir, "GSE165816_omnipath_lr_component_expression_by_sample_celltype.tsv")
write.table(expr_summary, expr_out, sep = "\t", quote = FALSE, row.names = FALSE)

avg_mat <- matrix(0, nrow = length(components_present), ncol = nrow(groups), dimnames = list(components_present, groups$group_id))
pct_mat <- matrix(0, nrow = length(components_present), ncol = nrow(groups), dimnames = list(components_present, groups$group_id))
avg_mat[cbind(expr_summary$gene, as.character(expr_summary$group_id))] <- expr_summary$avg_expr
pct_mat[cbind(expr_summary$gene, as.character(expr_summary$group_id))] <- expr_summary$pct_expr

score_side <- function(components_vec, mat_values) {
  out <- matrix(0, nrow = length(components_vec), ncol = ncol(mat_values))
  for (i in seq_along(components_vec)) {
    g <- strsplit(components_vec[[i]], ";", fixed = TRUE)[[1]]
    if (length(g) == 1) {
      out[i, ] <- mat_values[g, ]
    } else {
      out[i, ] <- apply(mat_values[g, , drop = FALSE], 2, min)
    }
  }
  out
}

message("Scoring ligand and receptor sides...")
lig_avg <- score_side(lr$source_components, avg_mat)
lig_pct <- score_side(lr$source_components, pct_mat)
rec_avg <- score_side(lr$target_components, avg_mat)
rec_pct <- score_side(lr$target_components, pct_mat)

dimnames(lig_avg) <- list(lr$lr_index, groups$group_id)
dimnames(lig_pct) <- dimnames(lig_avg)
dimnames(rec_avg) <- dimnames(lig_avg)
dimnames(rec_pct) <- dimnames(lig_avg)

long_side <- function(avg, pct, side_name) {
  grid <- expand.grid(lr_index = lr$lr_index, group_id = groups$group_id)
  grid[[paste0(side_name, "_avg")]] <- as.vector(avg)
  grid[[paste0(side_name, "_pct")]] <- as.vector(pct)
  merge(grid, groups, by = "group_id", all.x = TRUE)
}

lig <- long_side(lig_avg, lig_pct, "ligand")
names(lig)[names(lig) == "broad_cell_type"] <- "sender"
lig <- lig[, c("lr_index", "sample_id", "geo_accession", "disease", "sender", "ligand_avg", "ligand_pct")]

rec <- long_side(rec_avg, rec_pct, "receptor")
names(rec)[names(rec) == "broad_cell_type"] <- "receiver"
rec <- rec[, c("lr_index", "sample_id", "disease", "receiver", "receptor_avg", "receptor_pct")]

message("Combining sender and receiver cell types...")
scores <- merge(lig, rec, by = c("lr_index", "sample_id", "disease"), allow.cartesian = TRUE)
scores$score <- scores$ligand_avg * scores$ligand_pct * scores$receptor_avg * scores$receptor_pct
scores$key <- paste(scores$lr_index, scores$sender, scores$receiver, sep = "\r")
scores <- merge(
  scores,
  lr[, c(
    "lr_index", "source_genesymbol", "target_genesymbol", "source_components",
    "target_components", "interaction", "sources", "references", "curation_effort"
  )],
  by = "lr_index",
  all.x = TRUE
)

score_out <- file.path(out_dir, "GSE165816_omnipath_lr_sample_scores.tsv")
write.table(scores, score_out, sep = "\t", quote = FALSE, row.names = FALSE)

message("Testing healer vs nonhealer interaction scores...")
mean_by_key <- aggregate(score ~ key + disease, scores, mean)
wide <- reshape(mean_by_key, idvar = "key", timevar = "disease", direction = "wide")
names(wide) <- sub("^score\\.", "mean_", names(wide))
wide[is.na(wide)] <- 0
wide$max_mean <- pmax(wide$`mean_DFU-healer`, wide$`mean_DFU-nonhealer`)
keys_to_test <- wide$key[wide$max_mean > 1e-8]
idx_by_key <- split(seq_len(nrow(scores)), scores$key)

test_list <- vector("list", length(keys_to_test))
for (i in seq_along(keys_to_test)) {
  k <- keys_to_test[[i]]
  z <- scores[idx_by_key[[k]], ]
  healer <- z$score[z$disease == "DFU-healer"]
  nonhealer <- z$score[z$disease == "DFU-nonhealer"]
  first <- z[1, ]
  test_list[[i]] <- data.frame(
    lr_index = first$lr_index,
    interaction = first$interaction,
    ligand = first$source_genesymbol,
    receptor = first$target_genesymbol,
    ligand_components = first$source_components,
    receptor_components = first$target_components,
    sender = first$sender,
    receiver = first$receiver,
    curation_effort = first$curation_effort,
    sources = first$sources,
    n_healer = length(healer),
    n_nonhealer = length(nonhealer),
    mean_healer = mean(healer),
    mean_nonhealer = mean(nonhealer),
    delta_healer_minus_nonhealer = mean(healer) - mean(nonhealer),
    ratio_healer_over_nonhealer = (mean(healer) + 1e-8) / (mean(nonhealer) + 1e-8),
    p_wilcox = tryCatch(wilcox.test(healer, nonhealer, exact = FALSE)$p.value, error = function(e) NA_real_)
  )
}
tests <- do.call(rbind, test_list)
tests$padj_bh <- p.adjust(tests$p_wilcox, method = "BH")
tests <- tests[order(tests$p_wilcox, -tests$delta_healer_minus_nonhealer), ]

test_out <- file.path(out_dir, "GSE165816_omnipath_lr_healer_vs_nonhealer_tests.tsv")
write.table(tests, test_out, sep = "\t", quote = FALSE, row.names = FALSE)

focus_sender <- c("fibroblast_stromal", "myeloid", "endothelial")
focus_receiver <- c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte")
top_healer <- tests[
  tests$delta_healer_minus_nonhealer > 0 &
    tests$sender %in% focus_sender &
    tests$receiver %in% focus_receiver,
]
top_healer <- top_healer[order(top_healer$p_wilcox, -top_healer$ratio_healer_over_nonhealer), ]
top_healer <- head(top_healer, 100)

top_out <- file.path(out_dir, "GSE165816_omnipath_lr_top100_healer_interactions.tsv")
write.table(top_healer, top_out, sep = "\t", quote = FALSE, row.names = FALSE)

axis_genes <- c("IL11", "IL6", "THBS1", "TNC", "SERPINE1", "INHBA", "WNT5A", "CXCL5", "CXCL8", "CCL20", "FPR2", "VEGFA", "PDGFA", "TNFSF12")
axis_subset <- tests[
  tests$ligand %in% axis_genes |
    tests$receptor %in% axis_genes |
    grepl(paste(axis_genes, collapse = "|"), tests$ligand_components) |
    grepl(paste(axis_genes, collapse = "|"), tests$receptor_components),
]
axis_subset$padj_candidate_bh <- p.adjust(axis_subset$p_wilcox, method = "BH")
axis_subset <- axis_subset[order(axis_subset$p_wilcox, -axis_subset$delta_healer_minus_nonhealer), ]
axis_out <- file.path(out_dir, "GSE165816_omnipath_lr_candidate_axis_tests.tsv")
write.table(axis_subset, axis_out, sep = "\t", quote = FALSE, row.names = FALSE)

summary_pair <- aggregate(
  cbind(
    healer_higher_raw_p05 = tests$delta_healer_minus_nonhealer > 0 & tests$p_wilcox < 0.05,
    nonhealer_higher_raw_p05 = tests$delta_healer_minus_nonhealer < 0 & tests$p_wilcox < 0.05
  ) ~ sender + receiver,
  tests,
  sum
)
summary_out <- file.path(out_dir, "GSE165816_omnipath_lr_sender_receiver_summary.tsv")
write.table(summary_pair, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

plot_top <- head(top_healer, 30)
if (nrow(plot_top) > 0) {
  plot_top$label <- paste(plot_top$sender, "->", plot_top$receiver, plot_top$interaction)
  plot_top$label <- factor(plot_top$label, levels = rev(plot_top$label))
  p <- ggplot(plot_top, aes(x = label, y = delta_healer_minus_nonhealer, fill = sender)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Mean LR score delta, healer - nonhealer", fill = "Sender") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
  ggsave(
    file.path(fig_dir, "GSE165816_omnipath_lr_top30_healer_interactions.png"),
    p,
    width = 9,
    height = 7,
    dpi = 220
  )
}

summary_long <- rbind(
  data.frame(sender = summary_pair$sender, receiver = summary_pair$receiver, direction = "healer_higher", n = summary_pair$healer_higher_raw_p05),
  data.frame(sender = summary_pair$sender, receiver = summary_pair$receiver, direction = "nonhealer_higher", n = summary_pair$nonhealer_higher_raw_p05)
)
p2 <- ggplot(summary_long, aes(x = sender, y = receiver, fill = n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~ direction) +
  scale_fill_gradient(low = "#F5F5F5", high = "#0072B2") +
  labs(x = "Sender", y = "Receiver", fill = "Raw p<0.05\ninteractions") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(
  file.path(fig_dir, "GSE165816_omnipath_lr_sender_receiver_heatmap.png"),
  p2,
  width = 8,
  height = 4.8,
  dpi = 220
)

message("Wrote: ", expr_out)
message("Wrote: ", score_out)
message("Wrote: ", test_out)
message("Wrote: ", top_out)
message("Wrote: ", axis_out)
message("Wrote: ", summary_out)
