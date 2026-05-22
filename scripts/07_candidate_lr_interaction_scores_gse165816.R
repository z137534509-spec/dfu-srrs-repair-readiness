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
out_dir <- file.path(project_dir, "results", "candidate_lr_gse165816")
fig_dir <- file.path(project_dir, "figures", "candidate_lr_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading Seurat object...")
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

get_data_matrix <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
  )
}

lr <- data.frame(
  axis = c(
    "THBS1_integrin_CD47", "THBS1_integrin_CD47", "THBS1_integrin_CD47", "THBS1_LRP1",
    "IL6_gp130", "IL6_gp130",
    "IL11_gp130", "IL11_gp130",
    "INHBA_activin", "INHBA_activin", "INHBA_activin",
    "WNT5A_noncanonical", "WNT5A_noncanonical", "WNT5A_noncanonical", "WNT5A_noncanonical",
    "CCL20_CCR6",
    "CXCL5_CXCR", "CXCL5_CXCR", "CXCL8_CXCR", "CXCL8_CXCR",
    "ANXA1_FPR2",
    "TNC_integrin_TLR", "TNC_integrin_TLR", "TNC_integrin_TLR", "TNC_integrin_TLR",
    "SERPINE1_LRP1_PLAUR", "SERPINE1_LRP1_PLAUR",
    "TNFSF12_TWEAK", "TGFB_TGFBR", "TGFB_TGFBR",
    "VEGFA_KDR_FLT1", "VEGFA_KDR_FLT1",
    "PDGFA_PDGFRA_B", "PDGFA_PDGFRA_B"
  ),
  ligand = c(
    "THBS1", "THBS1", "THBS1", "THBS1",
    "IL6", "IL6",
    "IL11", "IL11",
    "INHBA", "INHBA", "INHBA",
    "WNT5A", "WNT5A", "WNT5A", "WNT5A",
    "CCL20",
    "CXCL5", "CXCL5", "CXCL8", "CXCL8",
    "ANXA1",
    "TNC", "TNC", "TNC", "TNC",
    "SERPINE1", "SERPINE1",
    "TNFSF12", "TGFB1", "TGFB2",
    "VEGFA", "VEGFA",
    "PDGFA", "PDGFA"
  ),
  receptor = c(
    "CD47", "ITGAV", "ITGB1", "LRP1",
    "IL6R", "IL6ST",
    "IL11RA", "IL6ST",
    "ACVR2A", "ACVR2B", "ACVR1B",
    "FZD2", "FZD5", "ROR1", "ROR2",
    "CCR6",
    "CXCR1", "CXCR2", "CXCR1", "CXCR2",
    "FPR2",
    "ITGAV", "ITGB1", "TLR4", "SDC4",
    "LRP1", "PLAUR",
    "TNFRSF12A", "TGFBR1", "TGFBR2",
    "KDR", "FLT1",
    "PDGFRA", "PDGFRB"
  ),
  stringsAsFactors = FALSE
)
lr$interaction <- paste(lr$ligand, lr$receptor, sep = "_")

focus_celltypes <- c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte", "pericyte_smc")
focus_diseases <- c("DFU-healer", "DFU-nonhealer")

genes <- unique(c(lr$ligand, lr$receptor))
mat <- get_data_matrix(obj)
genes_present <- intersect(genes, rownames(mat))
missing_genes <- setdiff(genes, genes_present)
message("Candidate LR genes present: ", length(genes_present), "/", length(genes))
if (length(missing_genes) > 0) {
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
}

meta <- obj@meta.data
meta$cell_id <- rownames(meta)
meta <- meta[meta$disease %in% focus_diseases & meta$broad_cell_type %in% focus_celltypes, ]
meta$sample_id <- if ("sample_code" %in% names(meta)) meta$sample_code else meta$geo_accession

groups <- unique(meta[, c("sample_id", "geo_accession", "disease", "broad_cell_type")])
groups <- groups[order(groups$disease, groups$sample_id, groups$broad_cell_type), ]

message("Summarizing expression by sample and cell type...")
expr_summary_list <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  g <- groups[i, ]
  cells <- meta$cell_id[
    meta$sample_id == g$sample_id &
      meta$disease == g$disease &
      meta$broad_cell_type == g$broad_cell_type
  ]
  m <- mat[genes_present, cells, drop = FALSE]
  expr_summary_list[[i]] <- data.frame(
    sample_id = g$sample_id,
    geo_accession = g$geo_accession,
    disease = g$disease,
    broad_cell_type = g$broad_cell_type,
    gene = genes_present,
    avg_expr = as.numeric(Matrix::rowMeans(m)),
    pct_expr = as.numeric(Matrix::rowMeans(m > 0)),
    cells = length(cells)
  )
}
expr_summary <- do.call(rbind, expr_summary_list)

expr_out <- file.path(out_dir, "GSE165816_candidate_lr_gene_expression_by_sample_celltype.tsv")
write.table(expr_summary, expr_out, sep = "\t", quote = FALSE, row.names = FALSE)

samples <- unique(meta[, c("sample_id", "geo_accession", "disease")])
sample_cell_grid <- expand.grid(
  sample_id = samples$sample_id,
  sender = focus_celltypes,
  receiver = focus_celltypes,
  lr_index = seq_len(nrow(lr)),
  stringsAsFactors = FALSE
)
sample_cell_grid <- merge(sample_cell_grid, samples, by = "sample_id", all.x = TRUE)

make_key <- function(sample_id, celltype, gene) paste(sample_id, celltype, gene, sep = "\r")
avg_lookup <- setNames(expr_summary$avg_expr, make_key(expr_summary$sample_id, expr_summary$broad_cell_type, expr_summary$gene))
pct_lookup <- setNames(expr_summary$pct_expr, make_key(expr_summary$sample_id, expr_summary$broad_cell_type, expr_summary$gene))
cell_lookup <- setNames(expr_summary$cells, make_key(expr_summary$sample_id, expr_summary$broad_cell_type, expr_summary$gene))

message("Scoring candidate interactions...")
scores <- sample_cell_grid
scores$axis <- lr$axis[scores$lr_index]
scores$ligand <- lr$ligand[scores$lr_index]
scores$receptor <- lr$receptor[scores$lr_index]
scores$interaction <- lr$interaction[scores$lr_index]
scores$ligand_avg <- avg_lookup[make_key(scores$sample_id, scores$sender, scores$ligand)]
scores$ligand_pct <- pct_lookup[make_key(scores$sample_id, scores$sender, scores$ligand)]
scores$receptor_avg <- avg_lookup[make_key(scores$sample_id, scores$receiver, scores$receptor)]
scores$receptor_pct <- pct_lookup[make_key(scores$sample_id, scores$receiver, scores$receptor)]
scores$sender_cells <- cell_lookup[make_key(scores$sample_id, scores$sender, scores$ligand)]
scores$receiver_cells <- cell_lookup[make_key(scores$sample_id, scores$receiver, scores$receptor)]

for (nm in c("ligand_avg", "ligand_pct", "receptor_avg", "receptor_pct", "sender_cells", "receiver_cells")) {
  scores[[nm]][is.na(scores[[nm]])] <- 0
}
scores$score <- scores$ligand_avg * scores$ligand_pct * scores$receptor_avg * scores$receptor_pct
scores <- scores[scores$ligand %in% genes_present & scores$receptor %in% genes_present, ]

score_out <- file.path(out_dir, "GSE165816_candidate_lr_sample_scores.tsv")
write.table(scores, score_out, sep = "\t", quote = FALSE, row.names = FALSE)

interaction_keys <- unique(scores[, c("axis", "interaction", "ligand", "receptor", "sender", "receiver")])
test_one <- function(i) {
  k <- interaction_keys[i, ]
  z <- scores[
    scores$axis == k$axis &
      scores$interaction == k$interaction &
      scores$sender == k$sender &
      scores$receiver == k$receiver,
  ]
  healer <- z$score[z$disease == "DFU-healer"]
  nonhealer <- z$score[z$disease == "DFU-nonhealer"]
  data.frame(
    axis = k$axis,
    interaction = k$interaction,
    ligand = k$ligand,
    receptor = k$receptor,
    sender = k$sender,
    receiver = k$receiver,
    n_healer = length(healer),
    n_nonhealer = length(nonhealer),
    mean_healer = mean(healer),
    mean_nonhealer = mean(nonhealer),
    delta_healer_minus_nonhealer = mean(healer) - mean(nonhealer),
    ratio_healer_over_nonhealer = (mean(healer) + 1e-6) / (mean(nonhealer) + 1e-6),
    p_wilcox = tryCatch(wilcox.test(healer, nonhealer, exact = FALSE)$p.value, error = function(e) NA_real_)
  )
}

tests <- do.call(rbind, lapply(seq_len(nrow(interaction_keys)), test_one))
tests$padj_bh <- p.adjust(tests$p_wilcox, method = "BH")
tests <- tests[order(tests$p_wilcox, -tests$delta_healer_minus_nonhealer), ]

test_out <- file.path(out_dir, "GSE165816_candidate_lr_healer_vs_nonhealer_tests.tsv")
write.table(tests, test_out, sep = "\t", quote = FALSE, row.names = FALSE)

focus_tests <- tests[
  tests$delta_healer_minus_nonhealer > 0 &
    tests$mean_healer > 0 &
    tests$sender %in% c("fibroblast_stromal", "myeloid", "endothelial") &
    tests$receiver %in% c("fibroblast_stromal", "myeloid", "endothelial", "keratinocyte"),
]
focus_tests <- head(focus_tests[order(focus_tests$p_wilcox, -focus_tests$ratio_healer_over_nonhealer), ], 30)
focus_out <- file.path(out_dir, "GSE165816_candidate_lr_top_healer_interactions.tsv")
write.table(focus_tests, focus_out, sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(focus_tests) > 0) {
  focus_tests$label <- paste(focus_tests$sender, "->", focus_tests$receiver, focus_tests$interaction)
  focus_tests$label <- factor(focus_tests$label, levels = rev(focus_tests$label))
  p <- ggplot(focus_tests, aes(x = label, y = delta_healer_minus_nonhealer, fill = axis)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(x = NULL, y = "Mean interaction score delta, healer - nonhealer", fill = "Axis") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
  ggsave(
    file.path(fig_dir, "GSE165816_candidate_lr_top_healer_interactions.png"),
    p,
    width = 9,
    height = 6.5,
    dpi = 220
  )
}

message("Wrote: ", expr_out)
message("Wrote: ", score_out)
message("Wrote: ", test_out)
message("Wrote: ", focus_out)
