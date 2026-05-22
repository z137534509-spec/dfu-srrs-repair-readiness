options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(igraph)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
obj_file <- file.path(project_dir, "results", "seurat_gse165816_foot_skin", "GSE165816_foot_skin_seurat_firstpass.rds")
out_dir <- file.path(project_dir, "results", "fibroblast_trajectory_gse165816")
fig_dir <- file.path(project_dir, "figures", "fibroblast_trajectory_gse165816")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"

fib <- subset(obj, subset = broad_cell_type == "fibroblast_stromal")
fib <- NormalizeData(fib, verbose = FALSE)
fib <- FindVariableFeatures(fib, nfeatures = 3000, verbose = FALSE)
fib <- ScaleData(fib, features = VariableFeatures(fib), verbose = FALSE)
fib <- RunPCA(fib, features = VariableFeatures(fib), npcs = 30, verbose = FALSE)
fib <- FindNeighbors(fib, dims = 1:20, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.55, verbose = FALSE)
fib <- RunUMAP(fib, dims = 1:20, verbose = FALSE)

program_sets <- list(
  fibroblast_sender_ligand_shortlist = c("IL11", "CCL20", "INHBA", "SERPINE1", "IL6", "THBS1", "TNC", "WNT5A", "ADAM12", "PTGS2"),
  fibroblast_repair_activation = c("THBS1", "WISP1", "IL11", "INHBA", "HIF1A", "TNFRSF12A", "CCL20", "MMP10", "COL10A1", "MMP3", "SERPINE1", "IL6", "PTGS2", "WNT5A", "TNFAIP6", "PDPN", "ITGA5", "ADAM12", "PLOD2"),
  published_healing_enriched_fibroblast = c("MMP1", "MMP3", "MMP11", "HIF1A", "CHI3L1", "TNFAIP6"),
  ecm_remodeling_migration = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1", "COL6A2", "COL10A1", "THBS1", "TNC", "MMP1", "MMP3", "MMP10", "SERPINE1", "ITGA5", "ITGAV", "PDPN", "TNFAIP6", "CTHRC1"),
  repair_inflammatory_signaling = c("IL6", "IL11", "INHBA", "AIM2", "CCL20", "CXCL5", "CXCL8", "TNFAIP3", "TNFRSF12A", "PTGS2", "FPR2", "TNIP3", "SLC39A8"),
  gp130_stromal = c("IL6", "IL11", "IL6ST", "STAT3", "SOCS3", "JUNB"),
  resting_matrix_fibroblast = c("DCN", "LUM", "CFD", "FBLN1", "COL1A1", "COL1A2", "COL3A1")
)

score_one <- function(object, genes, name) {
  genes <- intersect(genes, rownames(object))
  if (length(genes) < 2) {
    object[[name]] <- NA_real_
    return(object)
  }
  AddModuleScore(object, features = list(genes), name = paste0(name, "_tmp"), seed = 1, search = FALSE)
}

for (nm in names(program_sets)) {
  fib <- score_one(fib, program_sets[[nm]], nm)
  tmp_col <- paste0(nm, "_tmp1")
  if (tmp_col %in% colnames(fib@meta.data)) {
    fib[[nm]] <- fib@meta.data[[tmp_col]]
    fib@meta.data[[tmp_col]] <- NULL
  }
}

z <- function(x) as.numeric(scale(x))
fib$repair_axis <- rowMeans(
  cbind(
    z(fib$fibroblast_repair_activation),
    z(fib$published_healing_enriched_fibroblast),
    z(fib$fibroblast_sender_ligand_shortlist),
    z(fib$ecm_remodeling_migration),
    z(fib$repair_inflammatory_signaling),
    z(fib$gp130_stromal)
  ),
  na.rm = TRUE
)

cluster_col <- "seurat_clusters"
fib@meta.data$fibro_cluster <- as.character(fib@meta.data[[cluster_col]])

emb_pca <- Embeddings(fib, "pca")[, 1:15, drop = FALSE]
emb_umap <- Embeddings(fib, "umap")
cluster_ids <- sort(unique(fib$fibro_cluster))

cluster_centroids <- bind_rows(lapply(cluster_ids, function(cl) {
  cells <- WhichCells(fib, expression = fibro_cluster == cl)
  data.frame(
    fibro_cluster = cl,
    cell_n = length(cells),
    t(colMeans(emb_pca[cells, , drop = FALSE])),
    umap_1 = mean(emb_umap[cells, 1]),
    umap_2 = mean(emb_umap[cells, 2]),
    repair_axis = mean(fib$repair_axis[cells], na.rm = TRUE),
    published_he_fibro = mean(fib$published_healing_enriched_fibroblast[cells], na.rm = TRUE),
    repair_activation = mean(fib$fibroblast_repair_activation[cells], na.rm = TRUE),
    ligand_score = mean(fib$fibroblast_sender_ligand_shortlist[cells], na.rm = TRUE),
    gp130_score = mean(fib$gp130_stromal[cells], na.rm = TRUE),
    ecm_score = mean(fib$ecm_remodeling_migration[cells], na.rm = TRUE)
  )
}))

pca_cols <- grep("^PC_", colnames(cluster_centroids), value = TRUE)
dist_mat <- as.matrix(dist(cluster_centroids[, pca_cols, drop = FALSE]))
rownames(dist_mat) <- cluster_centroids$fibro_cluster
colnames(dist_mat) <- cluster_centroids$fibro_cluster
edge_pairs <- t(combn(cluster_centroids$fibro_cluster, 2))
edge_df <- data.frame(
  from = edge_pairs[, 1],
  to = edge_pairs[, 2],
  weight = dist_mat[cbind(edge_pairs[, 1], edge_pairs[, 2])]
)
g_full <- graph_from_data_frame(
  edge_df,
  directed = FALSE,
  vertices = data.frame(name = cluster_centroids$fibro_cluster)
)
g_mst <- mst(g_full, weights = E(g_full)$weight)

root_cluster <- cluster_centroids %>%
  arrange(repair_axis, published_he_fibro, gp130_score) %>%
  slice(1) %>%
  pull(fibro_cluster)

cluster_dist <- distances(g_mst, v = root_cluster, to = V(g_mst), weights = E(g_mst)$weight)
cluster_pseudotime <- data.frame(
  fibro_cluster = colnames(cluster_dist),
  pseudotime_raw = as.numeric(cluster_dist[1, ])
)
pt_range <- max(cluster_pseudotime$pseudotime_raw) - min(cluster_pseudotime$pseudotime_raw)
cluster_pseudotime$fibro_pseudotime <- if (is.finite(pt_range) && pt_range > 0) {
  (cluster_pseudotime$pseudotime_raw - min(cluster_pseudotime$pseudotime_raw)) / pt_range
} else {
  rep(0, nrow(cluster_pseudotime))
}

pt_map <- setNames(cluster_pseudotime$fibro_pseudotime, cluster_pseudotime$fibro_cluster)
raw_pt_map <- setNames(cluster_pseudotime$pseudotime_raw, cluster_pseudotime$fibro_cluster)
cluster_lookup <- as.character(fib@meta.data[["fibro_cluster"]])
fib@meta.data$fibro_pseudotime <- as.numeric(pt_map[cluster_lookup])
fib@meta.data$pseudotime_raw <- as.numeric(raw_pt_map[cluster_lookup])

mst_edges <- as_data_frame(g_mst, what = "edges") %>%
  transmute(from = as.character(from), to = as.character(to), weight = weight) %>%
  left_join(cluster_centroids %>% select(fibro_cluster, x = umap_1, y = umap_2), by = c("from" = "fibro_cluster")) %>%
  left_join(cluster_centroids %>% select(fibro_cluster, xend = umap_1, yend = umap_2), by = c("to" = "fibro_cluster"))

cluster_summary <- fib@meta.data %>%
  group_by(fibro_cluster) %>%
  summarise(
    cell_n = n(),
    pseudotime = mean(fibro_pseudotime, na.rm = TRUE),
    repair_axis = mean(repair_axis, na.rm = TRUE),
    published_he_fibro = mean(published_healing_enriched_fibroblast, na.rm = TRUE),
    repair_activation = mean(fibroblast_repair_activation, na.rm = TRUE),
    ligand_score = mean(fibroblast_sender_ligand_shortlist, na.rm = TRUE),
    gp130_score = mean(gp130_stromal, na.rm = TRUE),
    ecm_score = mean(ecm_remodeling_migration, na.rm = TRUE),
    healer_fraction = mean(disease == "DFU-healer"),
    nonhealer_fraction = mean(disease == "DFU-nonhealer"),
    diabetic_fraction = mean(disease == "Non-DFU Diabetic"),
    healthy_fraction = mean(disease == "Non-diabetic"),
    top_disease = names(sort(table(disease), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  arrange(pseudotime)

dfu_meta <- fib@meta.data %>%
  filter(disease %in% c("DFU-healer", "DFU-nonhealer")) %>%
  mutate(healing_status = ifelse(disease == "DFU-healer", "Healer", "Non-healer"))

late_cutoff <- quantile(fib$fibro_pseudotime, 0.75, na.rm = TRUE)
high_repair_cutoff <- quantile(fib$repair_axis, 0.75, na.rm = TRUE)

sample_scores <- dfu_meta %>%
  group_by(sample_code, healing_status) %>%
  summarise(
    fibroblast_cell_n = n(),
    mean_pseudotime = mean(fibro_pseudotime, na.rm = TRUE),
    late_pseudotime_fraction = mean(fibro_pseudotime >= late_cutoff, na.rm = TRUE),
    mean_repair_axis = mean(repair_axis, na.rm = TRUE),
    high_repair_fraction = mean(repair_axis >= high_repair_cutoff, na.rm = TRUE),
    mean_he_fibro_score = mean(published_healing_enriched_fibroblast, na.rm = TRUE),
    mean_ligand_score = mean(fibroblast_sender_ligand_shortlist, na.rm = TRUE),
    mean_gp130_score = mean(gp130_stromal, na.rm = TRUE),
    mean_ecm_score = mean(ecm_remodeling_migration, na.rm = TRUE),
    .groups = "drop"
  )

test_cols <- c(
  "mean_pseudotime", "late_pseudotime_fraction", "mean_repair_axis",
  "high_repair_fraction", "mean_he_fibro_score", "mean_ligand_score",
  "mean_gp130_score", "mean_ecm_score"
)

sample_tests <- bind_rows(lapply(test_cols, function(feature) {
  dat <- sample_scores %>% select(healing_status, value = all_of(feature))
  data.frame(
    feature = feature,
    healer_n = sum(dat$healing_status == "Healer"),
    nonhealer_n = sum(dat$healing_status == "Non-healer"),
    healer_mean = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE),
    nonhealer_mean = mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    delta_healer_minus_nonhealer = mean(dat$value[dat$healing_status == "Healer"], na.rm = TRUE) -
      mean(dat$value[dat$healing_status == "Non-healer"], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(value ~ healing_status, data = dat)$p.value)
  )
})) %>%
  mutate(BH = p.adjust(p_wilcox, method = "BH"))

program_cols <- c(
  "fibroblast_repair_activation", "published_healing_enriched_fibroblast",
  "fibroblast_sender_ligand_shortlist", "ecm_remodeling_migration",
  "repair_inflammatory_signaling", "gp130_stromal", "repair_axis"
)
program_cor <- bind_rows(lapply(program_cols, function(feature) {
  data.frame(
    feature = feature,
    spearman_rho = suppressWarnings(cor(fib$fibro_pseudotime, fib@meta.data[[feature]], method = "spearman", use = "pairwise.complete.obs")),
    p_spearman = suppressWarnings(cor.test(fib$fibro_pseudotime, fib@meta.data[[feature]], method = "spearman")$p.value)
  )
})) %>%
  mutate(BH = p.adjust(p_spearman, method = "BH"))

candidate_genes <- c("IL6", "IL11", "CCL20", "INHBA", "TGFB1", "TNC", "PTGS2", "THBS1", "SERPINE1", "ADAM12", "WNT5A", "MMP1", "MMP3", "MMP11", "HIF1A", "TNFAIP6")
expr_mat <- LayerData(fib, assay = "RNA", layer = "data")
candidate_genes_present <- intersect(candidate_genes, rownames(expr_mat))
avg_expr <- t(as.matrix(expr_mat[candidate_genes_present, , drop = FALSE])) %>%
  as.data.frame()
gene_cor <- bind_rows(lapply(colnames(avg_expr), function(gene) {
  data.frame(
    gene = gene,
    spearman_rho = suppressWarnings(cor(fib$fibro_pseudotime, avg_expr[[gene]], method = "spearman", use = "pairwise.complete.obs")),
    p_spearman = suppressWarnings(cor.test(fib$fibro_pseudotime, avg_expr[[gene]], method = "spearman")$p.value)
  )
})) %>%
  mutate(BH = p.adjust(p_spearman, method = "BH")) %>%
  arrange(desc(spearman_rho))

saveRDS(fib, file.path(out_dir, "GSE165816_fibroblast_trajectory_seurat.rds"))
write.table(cluster_summary, file.path(out_dir, "GSE165816_fibroblast_trajectory_cluster_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mst_edges, file.path(out_dir, "GSE165816_fibroblast_trajectory_mst_edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_scores, file.path(out_dir, "GSE165816_fibroblast_trajectory_sample_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_tests, file.path(out_dir, "GSE165816_fibroblast_trajectory_healer_vs_nonhealer_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(program_cor, file.path(out_dir, "GSE165816_fibroblast_program_pseudotime_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_cor, file.path(out_dir, "GSE165816_fibroblast_candidate_gene_pseudotime_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

umap_df <- Embeddings(fib, "umap") %>%
  as.data.frame() %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  mutate(
    cell = rownames(.),
    disease = fib$disease,
    fibro_cluster = fib$fibro_cluster,
    fibro_pseudotime = fib$fibro_pseudotime,
    repair_axis = fib$repair_axis,
    he_fibro_score = fib$published_healing_enriched_fibroblast
  )

p_umap_cluster <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = fibro_cluster)) +
  geom_point(size = 0.15, alpha = 0.7) +
  geom_segment(data = mst_edges, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, color = "black", linewidth = 0.45) +
  geom_text_repel(data = cluster_summary %>% left_join(cluster_centroids %>% select(fibro_cluster, umap_1, umap_2), by = "fibro_cluster"),
                  aes(x = umap_1, y = umap_2, label = fibro_cluster), inherit.aes = FALSE, size = 3) +
  labs(color = "Cluster") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

p_umap_disease <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = disease)) +
  geom_point(size = 0.15, alpha = 0.65) +
  scale_color_manual(values = c("DFU-healer" = "#0072B2", "DFU-nonhealer" = "#D55E00", "Non-DFU Diabetic" = "#999999", "Non-diabetic" = "#009E73")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

p_umap_pseudotime <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = fibro_pseudotime)) +
  geom_point(size = 0.15, alpha = 0.75) +
  scale_color_viridis_c(option = "plasma") +
  labs(color = "Pseudotime") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

p_umap_repair <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = repair_axis)) +
  geom_point(size = 0.15, alpha = 0.75) +
  scale_color_gradient2(low = "#D55E00", mid = "grey95", high = "#0072B2") +
  labs(color = "Repair axis") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "GSE165816_fibroblast_trajectory_umap_overview.png"),
       (p_umap_cluster | p_umap_disease) / (p_umap_pseudotime | p_umap_repair),
       width = 10, height = 8, dpi = 240)

p_sample <- sample_scores %>%
  pivot_longer(cols = all_of(test_cols), names_to = "feature", values_to = "value") %>%
  ggplot(aes(x = healing_status, y = value, color = healing_status)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(size = fibroblast_cell_n), position = position_jitter(width = 0.08, height = 0), alpha = 0.9) +
  facet_wrap(~feature, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Healer" = "#0072B2", "Non-healer" = "#D55E00")) +
  labs(x = NULL, y = "Sample-level score", color = NULL, size = "Fibroblast cells") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "GSE165816_fibroblast_trajectory_sample_level_healer_vs_nonhealer.png"),
       p_sample, width = 11, height = 6.5, dpi = 240)

trend_df <- fib@meta.data %>%
  select(fibro_pseudotime, all_of(program_cols)) %>%
  pivot_longer(-fibro_pseudotime, names_to = "feature", values_to = "value")
p_trend <- ggplot(trend_df, aes(x = fibro_pseudotime, y = value)) +
  geom_point(size = 0.08, alpha = 0.08, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#0072B2", linewidth = 0.8) +
  facet_wrap(~feature, scales = "free_y", ncol = 3) +
  labs(x = "Fibroblast pseudotime", y = "Program score") +
  theme_classic(base_size = 10)
ggsave(file.path(fig_dir, "GSE165816_fibroblast_programs_along_pseudotime.png"),
       p_trend, width = 9.5, height = 7, dpi = 240)

gene_expr_long <- t(as.matrix(expr_mat[candidate_genes_present, , drop = FALSE])) %>%
  as.data.frame() %>%
  mutate(fibro_pseudotime = fib$fibro_pseudotime) %>%
  pivot_longer(-fibro_pseudotime, names_to = "gene", values_to = "expression")
p_gene <- ggplot(gene_expr_long, aes(x = fibro_pseudotime, y = expression)) +
  geom_point(size = 0.08, alpha = 0.08, color = "grey40") +
  geom_smooth(method = "loess", se = FALSE, color = "#0072B2", linewidth = 0.75) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  labs(x = "Fibroblast pseudotime", y = "Log-normalized expression") +
  theme_classic(base_size = 10)
ggsave(file.path(fig_dir, "GSE165816_fibroblast_candidate_genes_along_pseudotime.png"),
       p_gene, width = 11, height = 8, dpi = 240)

message("Wrote fibroblast trajectory outputs.")
