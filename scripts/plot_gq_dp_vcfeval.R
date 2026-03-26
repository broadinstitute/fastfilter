#!/usr/bin/env Rscript
# Plots from metrics.tsv produced by gq_dp_vcfeval_scan.sh:
# - GQ-only scan: line plots of precision / recall / F1 vs GQ threshold
# - Legacy GQ×DP: heatmaps (if dp_thr column has variation)
# With a `sample` column, writes one set of PNGs per sample under <out_dir>/<sample>/.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript plot_gq_dp_vcfeval.R <metrics.tsv> [out_plots_dir]")
}

tsv <- args[[1]]
out_dir <- if (length(args) >= 2) args[[2]] else dirname(tsv)

if (!file.exists(tsv)) stop("file not found: ", tsv)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2: install.packages('ggplot2')")
}

theme_ff <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.35),
      panel.grid.major = ggplot2::element_line(color = "black", linewidth = 0.1),
      panel.grid.minor = ggplot2::element_line(color = "black", linewidth = 0.05),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(color = "black"),
      plot.subtitle = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(color = "black", face = "bold")
    )
}

raw <- read.table(tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
raw <- raw[!is.na(raw$gq_thr) & !grepl("^#", as.character(raw$gq_thr)), , drop = FALSE]

if (!"sample" %in% names(raw)) {
  raw$sample <- "sample"
}
if (!"dp_thr" %in% names(raw)) {
  raw$dp_thr <- NA_real_
}

num <- function(x) suppressWarnings(as.numeric(x))
raw$gq_thr <- num(raw$gq_thr)
raw$dp_thr <- num(raw$dp_thr)
raw$precision <- num(raw$precision)
raw$recall <- num(raw$recall)
raw$f1 <- num(raw$f1)

plot_heatmap <- function(df, col, title, out_path) {
  ok <- !is.na(df[[col]])
  if (!any(ok)) {
    warning("no finite values for ", col)
    return(invisible(NULL))
  }
  p <- ggplot2::ggplot(df[ok, ], ggplot2::aes(x = factor(dp_thr), y = factor(gq_thr), fill = .data[[col]])) +
    ggplot2::geom_tile(color = "black", linewidth = 0.12) +
    ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    ggplot2::labs(
      title = title,
      x = "DP threshold",
      y = "GQ threshold",
      fill = col
    ) +
    theme_ff() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  ggplot2::ggsave(out_path, p, width = 8, height = 6, dpi = 150, bg = "white")
  message("wrote ", out_path)
}

plot_lines_vs_gq <- function(df, title_suffix, pref) {
  df <- df[order(df$gq_thr), , drop = FALSE]
  m <- rbind(
    data.frame(gq_thr = df$gq_thr, y = df$precision, metric = "Precision"),
    data.frame(gq_thr = df$gq_thr, y = df$recall, metric = "Recall"),
    data.frame(gq_thr = df$gq_thr, y = df$f1, metric = "F1")
  )
  m <- m[is.finite(m$y), , drop = FALSE]
  if (nrow(m) == 0) {
    warning("no finite overall metrics for ", title_suffix)
    return(invisible(NULL))
  }
  p <- ggplot2::ggplot(m, ggplot2::aes(x = gq_thr, y = y, color = metric, group = metric)) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::scale_color_manual(values = c(Precision = "#1b9e77", Recall = "#d95f02", F1 = "#7570b3")) +
    ggplot2::labs(
      title = paste0("vcfeval overall metrics vs GQ threshold (", title_suffix, ")"),
      subtitle = "GQ-only filtering (no DP cutoff)",
      x = "Minimum GQ (fastfilter)",
      y = "Value"
    ) +
    theme_ff() +
    ggplot2::theme(legend.position = "none")
  out_path <- paste0(pref, "_metrics_vs_gq.png")
  ggplot2::ggsave(out_path, p, width = 10, height = 3.8, dpi = 150, bg = "white")
  message("wrote ", out_path)
}

samples <- unique(as.character(raw$sample))
for (s in samples) {
  sub <- raw[raw$sample == s, , drop = FALSE]
  pref <- file.path(out_dir, s, "vcfeval_grid")
  dir.create(dirname(pref), recursive = TRUE, showWarnings = FALSE)

  dp_ok <- sub$dp_thr[!is.na(sub$dp_thr)]
  use_heatmap <- length(dp_ok) > 0 && length(unique(dp_ok)) > 1

  if (use_heatmap) {
    plot_heatmap(sub, "f1", paste0("vcfeval F1 (", s, ")"), paste0(pref, "_f1.png"))
    plot_heatmap(sub, "precision", paste0("vcfeval Precision (", s, ")"), paste0(pref, "_precision.png"))
    plot_heatmap(sub, "recall", paste0("vcfeval Recall / sensitivity (", s, ")"), paste0(pref, "_recall.png"))
  } else {
    plot_lines_vs_gq(sub, s, pref)
  }
}
