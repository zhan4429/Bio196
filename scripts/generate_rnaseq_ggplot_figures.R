#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

set.seed(20260412)

outdir <- "images"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(filename, plot, width, height) {
  png_path <- file.path(outdir, paste0(filename, ".png"))
  pdf_path <- file.path(outdir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
}

base_theme <- theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 19, color = "#1D2733"),
    plot.subtitle = element_text(size = 14, color = "#52616B"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    axis.title = element_text(color = "#1D2733"),
    axis.text = element_text(color = "#1D2733"),
    legend.position = "bottom",
    legend.title = element_text(color = "#1D2733"),
    legend.text = element_text(color = "#1D2733")
  )

sample_info <- tibble(
  sample = c("CTRL_1", "CTRL_2", "CTRL_3", "TRT_1", "TRT_2", "TRT_3"),
  condition = rep(c("Control", "Treatment"), each = 3),
  batch = rep(c("Batch A", "Batch B", "Batch A"), 2)
)

n_genes <- 1200
gene_ids <- sprintf("gene_%04d", seq_len(n_genes))
base_mean <- rgamma(n_genes, shape = 2.2, rate = 0.05)
de_genes <- sample(seq_len(n_genes), 180)
true_lfc <- rep(0, n_genes)
true_lfc[de_genes] <- rnorm(length(de_genes), mean = 0, sd = 1.6)
true_lfc[de_genes[seq_len(90)]] <- abs(true_lfc[de_genes[seq_len(90)]]) + 0.8
true_lfc[de_genes[91:180]] <- -abs(true_lfc[de_genes[91:180]]) - 0.8

counts <- matrix(0, nrow = n_genes, ncol = nrow(sample_info))
for (j in seq_len(nrow(sample_info))) {
  condition_effect <- ifelse(sample_info$condition[j] == "Treatment", 2 ^ true_lfc, 1)
  batch_effect <- ifelse(sample_info$batch[j] == "Batch B", 1.12, 1)
  library_size <- runif(1, 0.88, 1.18)
  mu <- base_mean * condition_effect * batch_effect * library_size
  counts[, j] <- rnbinom(n_genes, mu = mu, size = 14)
}
rownames(counts) <- gene_ids
colnames(counts) <- sample_info$sample

log_counts <- log2(counts + 1)

# PCA plot --------------------------------------------------------------------
pca <- prcomp(t(log_counts), scale. = TRUE)
pca_df <- as_tibble(pca$x[, 1:2], rownames = "sample") |>
  left_join(sample_info, by = "sample")
pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 4.8, stroke = 1.1) +
  geom_text(aes(label = sample), nudge_y = 8, size = 4.4, show.legend = FALSE) +
  scale_color_manual(values = c(Control = "#2364AA", Treatment = "#D7263D")) +
  labs(
    title = "PCA of transformed RNA-seq counts",
    subtitle = "Samples should group by biology, not by hidden technical effects",
    x = paste0("PC1 (", pct[1], "% variance)"),
    y = paste0("PC2 (", pct[2], "% variance)"),
    color = "Condition",
    shape = "Batch"
  ) +
  base_theme

save_plot("09-rnaseq-ggplot-pca", p_pca, width = 8.4, height = 6.0)

# Sample distance heatmap -----------------------------------------------------
sample_dist <- as.matrix(dist(t(log_counts)))
hc <- hclust(as.dist(sample_dist))
sample_order <- hc$labels[hc$order]

dist_df <- as.data.frame(as.table(sample_dist)) |>
  as_tibble() |>
  rename(sample_a = Var1, sample_b = Var2, distance = Freq) |>
  mutate(
    sample_a = factor(sample_a, levels = sample_order),
    sample_b = factor(sample_b, levels = rev(sample_order))
  )

p_heatmap <- ggplot(dist_df, aes(sample_a, sample_b, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.7) +
  coord_fixed() +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B", name = "Distance") +
  labs(
    title = "Sample-to-sample distance heatmap",
    subtitle = "Replicates from the same condition should usually be more similar",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 19, color = "#1D2733"),
    plot.subtitle = element_text(size = 14, color = "#52616B"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "#1D2733"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

save_plot("09-rnaseq-ggplot-sample-distance", p_heatmap, width = 7.6, height = 6.4)

# Volcano plot ----------------------------------------------------------------
result_tbl <- tibble(
  gene_id = gene_ids,
  log2FoldChange = true_lfc + rnorm(n_genes, sd = 0.45),
  baseMean = rowMeans(counts)
) |>
  mutate(
    signal = abs(true_lfc) + log10(baseMean + 1) / 4,
    pvalue = pmin(1, 10 ^ (-(signal + rexp(n_genes, rate = 2.8)))),
    padj = p.adjust(pvalue, method = "BH"),
    status = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Higher in treatment",
      padj < 0.05 & log2FoldChange < -1 ~ "Lower in treatment",
      TRUE ~ "Not significant"
    ),
    neg_log10_p = -log10(pvalue)
  )

p_volcano <- ggplot(result_tbl, aes(log2FoldChange, neg_log10_p, color = status)) +
  geom_point(alpha = 0.72, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.45) +
  scale_color_manual(
    values = c(
      "Higher in treatment" = "#D7263D",
      "Lower in treatment" = "#2364AA",
      "Not significant" = "grey72"
    )
  ) +
  labs(
    title = "Volcano plot of differential expression results",
    subtitle = "Large fold changes with small p-values appear in the upper corners",
    x = "log2 fold change",
    y = "-log10 p-value",
    color = NULL
  ) +
  base_theme

save_plot("09-rnaseq-ggplot-volcano", p_volcano, width = 8.4, height = 6.0)

message("Generated RNA-seq ggplot figures in ", normalizePath(outdir))
