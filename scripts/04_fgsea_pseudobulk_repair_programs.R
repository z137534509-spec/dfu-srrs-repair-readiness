options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
de_dir <- file.path(project_dir, "results", "pseudobulk_de_gse165816")
out_dir <- file.path(project_dir, "results", "fgsea_pseudobulk_gse165816")
fig_dir <- file.path(project_dir, "figures", "fgsea_pseudobulk_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

message("Project: ", project_dir)

read_msig <- function(collection, subcollection = NULL) {
  dat <- tryCatch(
    {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(species = "Homo sapiens", category = collection)
      } else {
        msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
      }
    },
    error = function(e) {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
      } else {
        msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
      }
    }
  )
  split(dat$gene_symbol, dat$gs_name)
}

hallmark_sets <- read_msig("H")
go_bp_sets <- read_msig("C5", "GO:BP")

custom_sets <- list(
  CUSTOM_DFU_FIBROBLAST_REPAIR_ACTIVATION = c(
    "THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20",
    "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A",
    "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2"
  ),
  CUSTOM_ECM_REMODELING_MIGRATION = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
    "COL6A1", "COL6A2", "COL10A1", "THBS1", "TNC", "MMP1", "MMP3",
    "MMP10", "SERPINE1", "ITGA5", "ITGAV", "PDPN", "TNFAIP6", "CTHRC1"
  ),
  CUSTOM_REPAIR_INFLAMMATORY_SIGNALING = c(
    "IL6", "IL11", "INHBA", "AIM2", "CCL20", "CXCL5", "CXCL8",
    "TNFAIP3", "TNFRSF12A", "PTGS2", "FPR2", "TNIP3", "SLC39A8"
  ),
  CUSTOM_RESOLUTION_METABOLIC_MYELOID = c(
    "FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8",
    "CES1", "CXCL5", "TNC", "COL4A1", "COL4A2"
  )
)

pathways <- c(hallmark_sets, go_bp_sets, custom_sets)

de_file <- file.path(de_dir, "GSE165816_all_focus_celltypes_DESeq2_healer_vs_nonhealer.tsv")
de <- read.delim(de_file, check.names = FALSE)
de <- de[is.finite(de$stat) & !is.na(de$gene), ]

run_one <- function(cell_type) {
  sub <- de[de$cell_type == cell_type, ]
  sub <- sub[order(abs(sub$stat), decreasing = TRUE), ]
  sub <- sub[!duplicated(sub$gene), ]
  ranks <- sub$stat
  names(ranks) <- sub$gene
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 5,
    maxSize = 500
  )
  fg <- as.data.frame(fg)
  fg$cell_type <- cell_type
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg[order(fg$padj, -abs(fg$NES)), ]
}

cell_types <- sort(unique(de$cell_type))
all_fgsea <- do.call(rbind, lapply(cell_types, run_one))
all_fgsea$direction <- ifelse(all_fgsea$NES > 0, "higher_in_healer", "higher_in_nonhealer")

all_out <- file.path(out_dir, "GSE165816_fgsea_hallmark_go_custom_by_celltype.tsv")
write.table(all_fgsea, all_out, sep = "\t", quote = FALSE, row.names = FALSE)

focus_pattern <- paste(
  c(
    "CUSTOM_", "EXTRACELLULAR", "ECM", "COLLAGEN", "MATRIX", "WOUND",
    "MIGRATION", "MOTILITY", "ANGIOGENESIS", "VASCUL", "INFLAM",
    "CYTOKINE", "INTERLEUKIN", "CHEMOKINE", "TGF", "WNT", "HYPOX",
    "EPITHELIAL", "MESENCHYM", "RESPONSE_TO_WOUNDING"
  ),
  collapse = "|"
)
focus <- all_fgsea[grepl(focus_pattern, all_fgsea$pathway, ignore.case = TRUE), ]
focus <- focus[order(focus$cell_type, focus$padj, -abs(focus$NES)), ]
focus_out <- file.path(out_dir, "GSE165816_fgsea_focus_repair_terms.tsv")
write.table(focus, focus_out, sep = "\t", quote = FALSE, row.names = FALSE)

top_focus <- do.call(
  rbind,
  lapply(split(focus, focus$cell_type), function(x) {
    x <- x[x$padj <= 0.25, ]
    x <- x[order(x$direction, x$padj, -abs(x$NES)), ]
    head(x, 20)
  })
)
top_out <- file.path(out_dir, "GSE165816_fgsea_top_focus_terms_padj025.tsv")
write.table(top_focus, top_out, sep = "\t", quote = FALSE, row.names = FALSE)

plot_focus <- function(cell_type) {
  x <- focus[focus$cell_type == cell_type & focus$padj <= 0.25, ]
  if (nrow(x) == 0) {
    return(invisible(NULL))
  }
  x <- x[order(x$NES), ]
  x <- rbind(head(x, 12), tail(x, 12))
  x <- x[!duplicated(x$pathway), ]
  x$pathway_label <- gsub("^GOBP_|^HALLMARK_|^CUSTOM_", "", x$pathway)
  x$pathway_label <- gsub("_", " ", x$pathway_label)
  x$pathway_label <- factor(x$pathway_label, levels = x$pathway_label[order(x$NES)])

  p <- ggplot(x, aes(x = pathway_label, y = NES, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c(higher_in_healer = "#0072B2", higher_in_nonhealer = "#D55E00")) +
    labs(x = NULL, y = "Normalized enrichment score", fill = NULL, title = cell_type) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    filename = file.path(fig_dir, paste0("GSE165816_fgsea_focus_", cell_type, ".png")),
    plot = p,
    width = 8,
    height = max(4.5, 0.24 * nrow(x) + 1.5),
    dpi = 220
  )
}

invisible(lapply(cell_types, plot_focus))

summary_df <- data.frame(
  cell_type = cell_types,
  fgsea_terms = vapply(cell_types, function(ct) sum(all_fgsea$cell_type == ct), integer(1)),
  sig_padj_0_1_higher_healer = vapply(cell_types, function(ct) {
    sum(all_fgsea$cell_type == ct & all_fgsea$padj < 0.1 & all_fgsea$NES > 0, na.rm = TRUE)
  }, integer(1)),
  sig_padj_0_1_higher_nonhealer = vapply(cell_types, function(ct) {
    sum(all_fgsea$cell_type == ct & all_fgsea$padj < 0.1 & all_fgsea$NES < 0, na.rm = TRUE)
  }, integer(1)),
  focus_padj_0_25_higher_healer = vapply(cell_types, function(ct) {
    sum(focus$cell_type == ct & focus$padj < 0.25 & focus$NES > 0, na.rm = TRUE)
  }, integer(1)),
  focus_padj_0_25_higher_nonhealer = vapply(cell_types, function(ct) {
    sum(focus$cell_type == ct & focus$padj < 0.25 & focus$NES < 0, na.rm = TRUE)
  }, integer(1))
)
summary_out <- file.path(out_dir, "GSE165816_fgsea_summary_by_celltype.tsv")
write.table(summary_df, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote: ", all_out)
message("Wrote: ", focus_out)
message("Wrote: ", top_out)
message("Wrote: ", summary_out)
