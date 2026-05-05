#!/usr/bin/env Rscript
# Generate example figures for the ggplot2 chapter
# Uses simulated bioinformatics data

library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)

set.seed(42)
outdir <- "images"

# ── 1. Sample metadata ────────────────────────────────────────
samples <- tibble(
  sample_id = c("ctrl_1", "ctrl_2", "ctrl_3", "ctrl_4",
                "trt_1", "trt_2", "trt_3", "trt_4"),
  condition = factor(rep(c("control", "treatment"), each = 4),
                     levels = c("control", "treatment")),
  batch = rep(c("A", "B"), 4),
  reads_m = c(24.5, 27.1, 25.8, 26.3, 30.2, 33.5, 31.4, 29.8),
  mapped_percent = c(88.1, 90.4, 87.9, 89.2, 91.2, 89.5, 92.3, 90.8)
)

# ── 2. Simulated DE results ──────────────────────────────────
n_genes <- 15000
de <- tibble(
  gene_id = paste0("ENSG", sprintf("%05d", 1:n_genes)),
  symbol = c(
    # Known genes for top hits
    "TP53", "BRCA1", "MYC", "EGFR", "KRAS", "PIK3CA", "PTEN", "RB1",
    "APC", "VHL", "CDH1", "SMAD4", "NOTCH1", "FOXP3", "IL6",
    "TNF", "STAT3", "JAK2", "BCL2", "BRAF",
    paste0("Gene", 21:n_genes)
  ),
  baseMean = 10^runif(n_genes, 0.5, 4),
  log2FoldChange = rnorm(n_genes, 0, 0.8)
)

# Inject strong DE signal for some genes
up_idx <- 1:400
down_idx <- 401:750
de$log2FoldChange[up_idx] <- rnorm(length(up_idx), 2.5, 0.8)
de$log2FoldChange[down_idx] <- rnorm(length(down_idx), -2.2, 0.7)

# Generate p-values correlated with fold change
de <- de |>
  mutate(
    z_score = log2FoldChange / (0.3 + 0.5 / sqrt(baseMean)),
    pvalue = 2 * pnorm(-abs(z_score)),
    padj = p.adjust(pvalue, method = "BH"),
    direction = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    direction = factor(direction, levels = c("Up", "Down", "NS"))
  )

# ── 3. Simulated count data ─────────────────────────────────
gene_counts <- tibble(
  gene_id = paste0("gene", 1:8000),
  count = rnbinom(8000, mu = 80, size = 1.5)
)

# ── Fig 07c-scatter: Scatter plot ────────────────────────────
p_scatter <- ggplot(samples, aes(x = reads_m, y = mapped_percent,
                                  color = condition, shape = batch)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("control" = "#4575B4", "treatment" = "#D73027")) +
  theme_bw(base_size = 12) +
  labs(x = "Reads (millions)", y = "Mapped reads (%)",
       title = "Sequencing QC: Read Depth vs. Mapping Rate")

ggsave(file.path(outdir, "07c-ggplot2-scatter.png"), p_scatter,
       width = 6, height = 4, dpi = 200, bg = "white")

# ── Fig 07c-boxplot: Boxplot with jitter ─────────────────────
p_box <- ggplot(samples, aes(x = condition, y = reads_m, fill = condition)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2.5) +
  scale_fill_manual(values = c("control" = "#4575B4", "treatment" = "#D73027")) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Reads (millions)",
       title = "Read Depth by Condition") +
  theme(legend.position = "none")

ggsave(file.path(outdir, "07c-ggplot2-boxplot.png"), p_box,
       width = 4.5, height = 4, dpi = 200, bg = "white")

# ── Fig 07c-histogram: Gene count distribution ───────────────
p_hist <- ggplot(gene_counts, aes(x = count)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white", linewidth = 0.2) +
  scale_x_log10(labels = comma) +
  annotation_logticks(sides = "b") +
  theme_bw(base_size = 12) +
  labs(x = "Gene count (log scale)", y = "Number of genes",
       title = "Distribution of Gene Expression Counts")

ggsave(file.path(outdir, "07c-ggplot2-histogram.png"), p_hist,
       width = 6, height = 4, dpi = 200, bg = "white")

# ── Fig 07c-bar: DE gene counts ──────────────────────────────
de_summary <- de |>
  filter(direction != "NS") |>
  count(direction)

p_bar <- ggplot(de_summary, aes(x = direction, y = n, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Number of DE Genes",
       title = "Differentially Expressed Genes (padj < 0.05, |log2FC| > 1)") +
  theme(legend.position = "none")

ggsave(file.path(outdir, "07c-ggplot2-bar.png"), p_bar,
       width = 4.5, height = 4, dpi = 200, bg = "white")

# ── Fig 07c-volcano: Volcano plot ────────────────────────────
top_labels <- de |>
  filter(direction != "NS") |>
  slice_min(padj, n = 12)

p_volcano <- ggplot(de, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = direction), size = 0.4, alpha = 0.5) +
  scale_color_manual(
    values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70"),
    name = NULL
  ) +
  geom_text_repel(
    data = top_labels,
    aes(label = symbol),
    size = 3, max.overlaps = 20, segment.color = "grey50",
    min.segment.length = 0
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 12) +
  labs(
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value"),
    title = "Volcano Plot: Treatment vs. Control"
  ) +
  theme(legend.position = c(0.85, 0.85))

ggsave(file.path(outdir, "07c-ggplot2-volcano.png"), p_volcano,
       width = 7, height = 5, dpi = 200, bg = "white")

# ── Fig 07c-ma: MA plot ──────────────────────────────────────
p_ma <- ggplot(de, aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = direction), size = 0.4, alpha = 0.5) +
  scale_color_manual(
    values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70"),
    name = NULL
  ) +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw(base_size = 12) +
  labs(
    x = expression(log[10]~"Mean Expression"),
    y = expression(log[2]~"Fold Change"),
    title = "MA Plot"
  ) +
  theme(legend.position = c(0.85, 0.85))

ggsave(file.path(outdir, "07c-ggplot2-ma.png"), p_ma,
       width = 7, height = 5, dpi = 200, bg = "white")

# ── Fig 07c-pca: PCA plot (simulated) ────────────────────────
# Simulate PCA coordinates with condition separation
pca_data <- tibble(
  PC1 = c(rnorm(4, -15, 3), rnorm(4, 15, 3)),
  PC2 = c(rnorm(4, 2, 4), rnorm(4, -2, 4)),
  condition = samples$condition,
  batch = samples$batch,
  sample_id = samples$sample_id
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("control" = "#4575B4", "treatment" = "#D73027")) +
  theme_bw(base_size = 12) +
  labs(
    x = "PC1: 62% variance",
    y = "PC2: 14% variance",
    title = "PCA of Variance-Stabilized Counts"
  ) +
  coord_fixed(ratio = 0.5)

ggsave(file.path(outdir, "07c-ggplot2-pca.png"), p_pca,
       width = 6, height = 4.5, dpi = 200, bg = "white")

# ── Fig 07c-facet: Faceted histograms ────────────────────────
counts_per_sample <- tibble(
  sample_id = rep(samples$sample_id, each = 5000),
  count = unlist(lapply(1:8, function(i) {
    rnbinom(5000, mu = 60 + i * 10, size = 1.5)
  }))
)

p_facet <- ggplot(counts_per_sample |> filter(count > 0),
                  aes(x = log10(count))) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white", linewidth = 0.2) +
  facet_wrap(~sample_id, ncol = 4, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(x = expression(log[10]~"count"), y = "Genes",
       title = "Count Distribution per Sample")

ggsave(file.path(outdir, "07c-ggplot2-facet.png"), p_facet,
       width = 10, height = 5, dpi = 200, bg = "white")

# ── Fig 07c-patchwork: Multi-panel figure ────────────────────
p_a <- p_volcano +
  labs(title = NULL) +
  theme(legend.position = "none", plot.margin = margin(5, 10, 5, 5))

p_b <- p_ma +
  labs(title = NULL) +
  theme(legend.position = "none", plot.margin = margin(5, 5, 5, 10))

p_c <- de |>
  filter(direction != "NS") |>
  count(direction) |>
  ggplot(aes(x = direction, y = n, fill = direction)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "DE genes") +
  theme(legend.position = "none")

p_panel <- (p_a | p_b | p_c) +
  plot_layout(widths = c(3, 3, 1.5)) +
  plot_annotation(tag_levels = "A",
                  title = "Differential Expression: Treatment vs. Control",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(outdir, "07c-ggplot2-patchwork.png"), p_panel,
       width = 13, height = 4.5, dpi = 200, bg = "white")

# ── Fig 07c-theme: Theme comparison ─────────────────────────
p_base <- ggplot(samples, aes(x = condition, y = reads_m, fill = condition)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c("control" = "#4575B4", "treatment" = "#D73027")) +
  labs(x = NULL, y = "Reads (M)") +
  theme(legend.position = "none")

p_themes <- (p_base + theme_gray() + ggtitle("theme_gray()")) |
  (p_base + theme_bw() + ggtitle("theme_bw()")) |
  (p_base + theme_minimal() + ggtitle("theme_minimal()")) |
  (p_base + theme_classic() + ggtitle("theme_classic()"))

ggsave(file.path(outdir, "07c-ggplot2-themes.png"), p_themes,
       width = 12, height = 3.5, dpi = 200, bg = "white")

cat("All figures generated in", outdir, "\n")
