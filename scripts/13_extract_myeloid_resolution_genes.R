suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

my <- readRDS("results/myeloid_subcluster_gse165816/GSE165816_myeloid_subcluster_seurat.rds")
mat <- GetAssayData(my, assay = "RNA", layer = "data")
genes <- intersect(
  c(
    "FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CES1", "CXCL5",
    "TNC", "COL4A1", "COL4A2", "IL1B", "S100A8", "C1QA", "APOE", "CD74", "FCER1A"
  ),
  rownames(mat)
)
meta <- my@meta.data
out <- do.call(rbind, lapply(sort(unique(meta$myeloid_subcluster)), function(cl) {
  cells <- rownames(meta)[meta$myeloid_subcluster == cl]
  m <- mat[genes, cells, drop = FALSE]
  data.frame(
    myeloid_subcluster = cl,
    gene = genes,
    avg_expr = as.numeric(Matrix::rowMeans(m)),
    pct_expr = as.numeric(Matrix::rowMeans(m > 0))
  )
}))

write.table(
  out,
  "results/myeloid_subcluster_gse165816/GSE165816_myeloid_resolution_genes_by_subcluster.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(
  out %>%
    filter(gene %in% c("FPR2", "SLC7A11", "TNIP3", "STEAP4", "ACSL1", "SLC39A8", "CXCL5", "TNC")) %>%
    arrange(gene, desc(avg_expr)) %>%
    group_by(gene) %>%
    slice_head(n = 5)
)
