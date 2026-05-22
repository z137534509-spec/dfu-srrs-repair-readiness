options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
infile <- file.path(
  project_dir,
  "results", "seurat_gse165816_foot_skin",
  "GSE165816_foot_skin_sample_celltype_program_scores.tsv"
)
out_dir <- file.path(project_dir, "results", "composition_gse165816")
fig_dir <- file.path(project_dir, "figures", "composition_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(infile, check.names = FALSE)
x <- x[, c("geo_accession", "sample_code", "disease", "broad_cell_type", "cells")]

sample_meta <- unique(x[, c("geo_accession", "sample_code", "disease")])
cell_types <- sort(unique(x$broad_cell_type))
grid <- merge(sample_meta, data.frame(broad_cell_type = cell_types), all = TRUE)
merged <- merge(
  grid,
  x,
  by = c("geo_accession", "sample_code", "disease", "broad_cell_type"),
  all.x = TRUE
)
merged$cells[is.na(merged$cells)] <- 0L

totals <- aggregate(cells ~ geo_accession + sample_code + disease, merged, sum)
names(totals)[names(totals) == "cells"] <- "total_cells"
merged <- merge(merged, totals, by = c("geo_accession", "sample_code", "disease"))
merged$fraction <- merged$cells / merged$total_cells

prop_out <- file.path(out_dir, "GSE165816_sample_broad_celltype_fractions.tsv")
write.table(merged, prop_out, sep = "\t", quote = FALSE, row.names = FALSE)

dfu <- merged[merged$disease %in% c("DFU-healer", "DFU-nonhealer"), ]

test_one <- function(ct) {
  z <- dfu[dfu$broad_cell_type == ct, ]
  healer <- z$fraction[z$disease == "DFU-healer"]
  nonhealer <- z$fraction[z$disease == "DFU-nonhealer"]
  p <- tryCatch(wilcox.test(healer, nonhealer, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.frame(
    broad_cell_type = ct,
    n_healer = length(healer),
    n_nonhealer = length(nonhealer),
    mean_fraction_healer = mean(healer),
    mean_fraction_nonhealer = mean(nonhealer),
    median_fraction_healer = median(healer),
    median_fraction_nonhealer = median(nonhealer),
    delta_mean_healer_minus_nonhealer = mean(healer) - mean(nonhealer),
    p_wilcox = p
  )
}

tests <- do.call(rbind, lapply(cell_types, test_one))
tests$padj_bh <- p.adjust(tests$p_wilcox, method = "BH")
tests <- tests[order(tests$p_wilcox), ]

test_out <- file.path(out_dir, "GSE165816_DFU_healer_vs_nonhealer_celltype_fraction_wilcox.tsv")
write.table(tests, test_out, sep = "\t", quote = FALSE, row.names = FALSE)

focus_order <- c(
  "fibroblast_stromal", "keratinocyte", "myeloid", "endothelial",
  "pericyte_smc", "t_nk", "b_plasma", "schwann", "melanocyte"
)
merged$broad_cell_type <- factor(merged$broad_cell_type, levels = focus_order)
merged$sample_label <- paste(merged$sample_code, merged$disease, sep = " | ")
sample_order <- unique(merged[order(merged$disease, merged$sample_code), "sample_label"])
merged$sample_label <- factor(merged$sample_label, levels = sample_order)

palette <- c(
  fibroblast_stromal = "#0072B2",
  keratinocyte = "#D55E00",
  myeloid = "#009E73",
  endothelial = "#CC79A7",
  pericyte_smc = "#E69F00",
  t_nk = "#56B4E9",
  b_plasma = "#7F7F7F",
  schwann = "#8A63D2",
  melanocyte = "#111111"
)

p1 <- ggplot(merged, aes(x = sample_label, y = fraction, fill = broad_cell_type)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = palette, drop = FALSE) +
  coord_flip() +
  labs(x = NULL, y = "Fraction of QC-passed cells", fill = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
ggsave(
  file.path(fig_dir, "GSE165816_broad_celltype_fraction_by_sample.png"),
  p1,
  width = 9,
  height = 7,
  dpi = 220
)

dfu_focus <- dfu[dfu$broad_cell_type %in% focus_order[1:6], ]
dfu_focus$broad_cell_type <- factor(dfu_focus$broad_cell_type, levels = focus_order[1:6])
p2 <- ggplot(dfu_focus, aes(x = disease, y = fraction, fill = disease)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.8, alpha = 0.8) +
  facet_wrap(~ broad_cell_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(`DFU-healer` = "#0072B2", `DFU-nonhealer` = "#D55E00")) +
  labs(x = NULL, y = "Fraction of QC-passed cells", fill = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(
  file.path(fig_dir, "GSE165816_DFU_healer_vs_nonhealer_celltype_fraction_boxplot.png"),
  p2,
  width = 8,
  height = 5,
  dpi = 220
)

message("Wrote: ", prop_out)
message("Wrote: ", test_out)
