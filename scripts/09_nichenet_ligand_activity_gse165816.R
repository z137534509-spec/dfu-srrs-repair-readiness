options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(nichenetr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
seurat_file <- file.path(
  project_dir,
  "results", "seurat_gse165816_foot_skin",
  "GSE165816_foot_skin_seurat_firstpass.rds"
)
de_file <- file.path(
  project_dir,
  "results", "pseudobulk_de_gse165816",
  "GSE165816_all_focus_celltypes_DESeq2_healer_vs_nonhealer.tsv"
)
ltm_file <- file.path(project_dir, "data", "metadata", "nichenet_ligand_target_matrix_zenodo3260758.rds")
lr_file <- file.path(project_dir, "data", "metadata", "nichenet_lr_network_zenodo3260758.rds")

out_dir <- file.path(project_dir, "results", "nichenet_gse165816")
fig_dir <- file.path(project_dir, "figures", "nichenet_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

focus_celltypes <- c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte", "pericyte_smc")
focus_diseases <- c("DFU-healer", "DFU-nonhealer")
min_pct <- 0.05
min_avg <- 0.02
min_geneset <- 20
fallback_top_n <- 100

message("Loading NicheNet networks...")
ligand_target_matrix <- readRDS(ltm_file)
lr_network <- readRDS(lr_file)
lr_network <- lr_network[!duplicated(lr_network[, c("from", "to")]), c("from", "to", "source", "database")]

message("Loading Seurat object...")
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
mat <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)
meta <- meta[meta$disease %in% focus_diseases & meta$broad_cell_type %in% focus_celltypes, ]

message("Summarizing expressed genes by cell type...")
expr_by_celltype <- vector("list", length(focus_celltypes))
names(expr_by_celltype) <- focus_celltypes
for (ct in focus_celltypes) {
  cells <- meta$cell_id[meta$broad_cell_type == ct]
  m <- mat[, cells, drop = FALSE]
  expr_by_celltype[[ct]] <- data.frame(
    broad_cell_type = ct,
    gene = rownames(m),
    avg_expr = as.numeric(Matrix::rowMeans(m)),
    pct_expr = as.numeric(Matrix::rowMeans(m > 0)),
    cells = length(cells)
  )
}
expr_summary <- do.call(rbind, expr_by_celltype)
expr_out <- file.path(out_dir, "GSE165816_nichenet_expression_by_celltype.tsv")
write.table(expr_summary, expr_out, sep = "\t", quote = FALSE, row.names = FALSE)

expressed_genes <- lapply(expr_by_celltype, function(x) {
  x$gene[x$pct_expr >= min_pct & x$avg_expr >= min_avg]
})

de <- read.delim(de_file, check.names = FALSE)
de <- de[de$cell_type %in% focus_celltypes, ]
de$padj_filled <- ifelse(is.na(de$padj), 1, de$padj)

select_geneset <- function(receiver) {
  z <- de[de$cell_type == receiver & is.finite(de$stat), ]
  primary <- z$gene[z$log2FoldChange > 0 & z$padj_filled < 0.1]
  method <- "padj_lt_0.1_positive"
  if (length(primary) < min_geneset) {
    zpos <- z[z$log2FoldChange > 0, ]
    zpos <- zpos[order(zpos$pvalue, -zpos$stat), ]
    primary <- head(zpos$gene, fallback_top_n)
    method <- paste0("fallback_top_", fallback_top_n, "_positive_by_pvalue")
  }
  primary <- unique(primary)
  primary <- intersect(primary, rownames(ligand_target_matrix))
  background <- intersect(expressed_genes[[receiver]], rownames(ligand_target_matrix))
  list(geneset = primary, background = background, method = method)
}

run_pair <- function(sender, receiver) {
  selected <- select_geneset(receiver)
  geneset <- selected$geneset
  background <- selected$background
  sender_expressed <- expressed_genes[[sender]]
  receiver_expressed <- expressed_genes[[receiver]]
  potential_ligands <- unique(lr_network$from[
    lr_network$from %in% sender_expressed &
      lr_network$to %in% receiver_expressed
  ])
  potential_ligands <- intersect(potential_ligands, colnames(ligand_target_matrix))

  if (length(geneset) < 5 || length(background) < 50 || length(potential_ligands) < 5) {
    return(list(
      activities = data.frame(),
      links = data.frame(),
      meta = data.frame(
        sender = sender,
        receiver = receiver,
        geneset_n = length(geneset),
        background_n = length(background),
        potential_ligands_n = length(potential_ligands),
        geneset_method = selected$method,
        status = "skipped"
      )
    ))
  }

  activities <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  activities <- activities[order(-activities$pearson), ]
  activities$sender <- sender
  activities$receiver <- receiver
  activities$geneset_n <- length(geneset)
  activities$background_n <- length(background)
  activities$potential_ligands_n <- length(potential_ligands)
  activities$geneset_method <- selected$method

  top_ligands <- head(activities$test_ligand, 15)
  links <- do.call(rbind, lapply(top_ligands, function(lig) {
    x <- nichenetr::get_weighted_ligand_target_links(
      ligand = lig,
      geneset = geneset,
      ligand_target_matrix = ligand_target_matrix,
      n = 50
    )
    if (nrow(x) == 0) {
      return(data.frame())
    }
    x$sender <- sender
    x$receiver <- receiver
    x
  }))
  if (nrow(links) > 0) {
    links$geneset_method <- selected$method
  }

  list(
    activities = activities,
    links = links,
    meta = data.frame(
      sender = sender,
      receiver = receiver,
      geneset_n = length(geneset),
      background_n = length(background),
      potential_ligands_n = length(potential_ligands),
      geneset_method = selected$method,
      status = "completed"
    )
  )
}

pairs <- expand.grid(sender = focus_celltypes, receiver = focus_celltypes, stringsAsFactors = FALSE)
message("Running NicheNet ligand activity for ", nrow(pairs), " sender-receiver pairs...")
res <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  message("Pair ", i, "/", nrow(pairs), ": ", pairs$sender[i], " -> ", pairs$receiver[i])
  res[[i]] <- run_pair(pairs$sender[i], pairs$receiver[i])
}

activities <- do.call(rbind, lapply(res, `[[`, "activities"))
links <- do.call(rbind, lapply(res, `[[`, "links"))
run_meta <- do.call(rbind, lapply(res, `[[`, "meta"))

act_out <- file.path(out_dir, "GSE165816_nichenet_ligand_activities.tsv")
links_out <- file.path(out_dir, "GSE165816_nichenet_ligand_target_links.tsv")
meta_out <- file.path(out_dir, "GSE165816_nichenet_run_metadata.tsv")
write.table(activities, act_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(links, links_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(run_meta, meta_out, sep = "\t", quote = FALSE, row.names = FALSE)

priority_pairs <- data.frame(
  sender = c("fibroblast_stromal", "fibroblast_stromal", "fibroblast_stromal", "fibroblast_stromal", "myeloid", "endothelial", "pericyte_smc"),
  receiver = c("myeloid", "endothelial", "keratinocyte", "fibroblast_stromal", "fibroblast_stromal", "fibroblast_stromal", "fibroblast_stromal")
)
priority_key <- paste(priority_pairs$sender, priority_pairs$receiver, sep = "->")
activities$pair <- paste(activities$sender, activities$receiver, sep = "->")
priority <- activities[activities$pair %in% priority_key, ]
priority <- do.call(rbind, lapply(split(priority, priority$pair), function(x) head(x[order(-x$pearson), ], 20)))
priority_out <- file.path(out_dir, "GSE165816_nichenet_priority_pair_top_ligands.tsv")
write.table(priority, priority_out, sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(priority) > 0) {
  plot_df <- do.call(rbind, lapply(split(priority, priority$pair), function(x) head(x[order(-x$pearson), ], 8)))
  plot_df$label <- paste(plot_df$pair, plot_df$test_ligand, sep = " | ")
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label[order(plot_df$pearson)]))
  p <- ggplot(plot_df, aes(x = label, y = pearson, fill = pair)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "NicheNet ligand activity (Pearson)", fill = "Pair") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
  ggsave(
    file.path(fig_dir, "GSE165816_nichenet_priority_pair_top_ligands.png"),
    p,
    width = 9,
    height = 8,
    dpi = 220
  )
}

global_top <- do.call(rbind, lapply(split(activities, activities$pair), function(x) head(x[order(-x$pearson), ], 5)))
global_out <- file.path(out_dir, "GSE165816_nichenet_top5_ligands_by_pair.tsv")
write.table(global_top, global_out, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", expr_out)
message("Wrote: ", act_out)
message("Wrote: ", links_out)
message("Wrote: ", meta_out)
message("Wrote: ", priority_out)
message("Wrote: ", global_out)
