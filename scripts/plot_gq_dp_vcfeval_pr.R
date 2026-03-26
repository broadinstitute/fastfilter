#!/usr/bin/env Rscript
# Precision–recall from gq_dp_vcfeval_scan metrics.tsv (stratified SNV / insertion / deletion).
# Writes a combined plot with one translucent line per sample (faceted by variant class).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript plot_gq_dp_vcfeval_pr.R <metrics.tsv> [out_plots_dir]")
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
      panel.grid.major = ggplot2::element_line(color = "black", linewidth = 0.12),
      panel.grid.minor = ggplot2::element_line(color = "black", linewidth = 0.06),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(color = "black"),
      plot.subtitle = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(color = "black", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black", linewidth = 0.2)
    )
}

raw <- read.table(tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
req <- c(
  "gq_thr",
  "precision_snv", "recall_snv", "precision_ins", "recall_ins", "precision_del", "recall_del"
)
miss <- setdiff(req, names(raw))
if (length(miss) > 0) {
  stop("metrics.tsv missing columns (re-run scan script): ", paste(miss, collapse = ", "))
}

if (!"sample" %in% names(raw)) {
  raw$sample <- "sample"
}
if (!"dp_thr" %in% names(raw)) {
  raw$dp_thr <- NA_real_
}

num <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA"] <- NA_character_
  suppressWarnings(as.numeric(x))
}

raw$gq_thr <- num(raw$gq_thr)
raw$dp_thr <- num(raw$dp_thr)
for (nm in grep("^precision_|^recall_", names(raw), value = TRUE)) {
  raw[[nm]] <- num(raw[[nm]])
}

long <- rbind(
  data.frame(
    sample = raw$sample, gq_thr = raw$gq_thr, dp_thr = raw$dp_thr,
    precision = raw$precision_snv, recall = raw$recall_snv,
    variant_class = "SNV"
  ),
  data.frame(
    sample = raw$sample, gq_thr = raw$gq_thr, dp_thr = raw$dp_thr,
    precision = raw$precision_ins, recall = raw$recall_ins,
    variant_class = "Insertion"
  ),
  data.frame(
    sample = raw$sample, gq_thr = raw$gq_thr, dp_thr = raw$dp_thr,
    precision = raw$precision_del, recall = raw$recall_del,
    variant_class = "Deletion"
  )
)
long <- long[is.finite(long$precision) & is.finite(long$recall), , drop = FALSE]
if (nrow(long) == 0) stop("no finite precision/recall rows")

long$variant_class <- factor(long$variant_class, levels = c("SNV", "Insertion", "Deletion"))

dp_ok <- long$dp_thr[!is.na(long$dp_thr)]
gq_only <- length(dp_ok) == 0 || length(unique(dp_ok)) <= 1

if (gq_only) {
  long <- long[order(long$variant_class, long$sample, long$gq_thr), , drop = FALSE]
  p <- ggplot2::ggplot(long, ggplot2::aes(x = recall, y = precision, group = sample)) +
    ggplot2::geom_path(color = "black", linewidth = 0.55, alpha = 0.28, lineend = "round") +
    ggplot2::geom_point(ggplot2::aes(color = factor(gq_thr)), size = 1.8, alpha = 0.9) +
    ggplot2::facet_wrap(~variant_class, ncol = 3) +
    ggplot2::scale_color_viridis_d(option = "D", name = "GQ\nthreshold") +
    ggplot2::labs(
      title = "Precision vs recall (vcfeval)",
      subtitle = "",
      x = "Sensitivity",
      y = "Precision"
    ) +
    ggplot2::coord_cartesian(xlim = c(0.65, 1.00), ylim = c(0.75, 1.00)) +
    theme_ff()
} else {
  long <- long[order(long$variant_class, long$sample, long$dp_thr, long$gq_thr), , drop = FALSE]
  long$sample_dp <- interaction(long$sample, long$dp_thr, drop = TRUE)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = recall, y = precision, group = sample_dp)) +
    ggplot2::geom_path(color = "black", linewidth = 0.5, alpha = 0.22, lineend = "round") +
    ggplot2::facet_wrap(~variant_class, ncol = 3) +
    ggplot2::labs(
      title = "Precision vs recall (vcfeval), one line per sample×DP path",
      subtitle = "Translucent black lines show sample variability",
      x = "Recall (sensitivity)",
      y = "Precision"
    ) +
    ggplot2::coord_cartesian(xlim = c(0.65, 1.00), ylim = c(0.75, 1.00)) +
    theme_ff()
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "vcfeval_grid_pr_by_sample.png")
ggplot2::ggsave(out_path, p, width = 11, height = 4, dpi = 150, bg = "white")
message("wrote ", out_path)

# Additional SNV-only zoomed figure for fine-scale inspection.
snv <- long[long$variant_class == "SNV", , drop = FALSE]
if (nrow(snv) > 0) {
  if (gq_only) {
    snv <- snv[order(snv$sample, snv$gq_thr), , drop = FALSE]
    p_snv <- ggplot2::ggplot(snv, ggplot2::aes(x = recall, y = precision, group = sample)) +
      ggplot2::geom_path(color = "black", linewidth = 0.55, alpha = 0.28, lineend = "round") +
        ggplot2::geom_point(ggplot2::aes(color = factor(gq_thr)), size = 2.0, alpha = 0.95) +
        ggplot2::scale_color_viridis_d(option = "D", name = "GQ\nthreshold") +
      ggplot2::labs(
        title = "Precision vs recall (vcfeval)",
        subtitle = "",
        x = "Recall",
        y = "Precision"
      ) +
      ggplot2::coord_cartesian(xlim = c(0.95, 1.00), ylim = c(0.95, 1.00)) +
      theme_ff()
  } else {
    snv$sample_dp <- interaction(snv$sample, snv$dp_thr, drop = TRUE)
    snv <- snv[order(snv$sample, snv$dp_thr, snv$gq_thr), , drop = FALSE]
    p_snv <- ggplot2::ggplot(snv, ggplot2::aes(x = recall, y = precision, group = sample_dp)) +
      ggplot2::geom_path(color = "black", linewidth = 0.5, alpha = 0.22, lineend = "round") +
      ggplot2::labs(
        title = "Precision vs recall",
        subtitle = "",
        x = "Recall",
        y = "Precision"
      ) +
      ggplot2::coord_cartesian(xlim = c(0.95, 1.00), ylim = c(0.95, 1.00)) +
      theme_ff()
  }
  out_snv <- file.path(out_dir, "vcfeval_grid_pr_by_sample_snv_zoom.png")
  ggplot2::ggsave(out_snv, p_snv, width = 6.5, height = 5, dpi = 150, bg = "white")
  message("wrote ", out_snv)
}

# Base-R 3-panel faceted-style figure with per-panel limits.
# Matches "plot(..., bty='n')" style and allows independent SNV zoom.
draw_panel <- function(df, cls, xlim, ylim, main_title) {
  sub <- df[df$variant_class == cls, , drop = FALSE]
  plot(
    NA, NA,
    xlim = xlim, ylim = ylim,
    xlab = "Recall",
    ylab = "Precision",
    main = main_title,
    bty = "n", xaxs = "i", yaxs = "i"
  )
  if (nrow(sub) == 0) return(invisible(NULL))

  if (cls == "SNV" && gq_only) {
    legend(
      "topleft",
      legend = paste0("GQ=", gq_levels),
      col = pal2,
      pch = 16,
      cex = 0.7,
      bty = "n"
    )
  }

  samps <- unique(sub$sample)
  for (s in samps) {
    d <- sub[sub$sample == s, , drop = FALSE]
    if (gq_only) {
      d <- d[order(d$gq_thr), , drop = FALSE]
    } else {
      d <- d[order(d$dp_thr, d$gq_thr), , drop = FALSE]
    }
    lines(d$recall, d$precision, col = rgb(0, 0, 0, alpha = 0.28), lwd = 1)

    if (gq_only) {
      cols <- pal2[match(d$gq_thr, gq_levels)]
      points(d$recall, d$precision, pch = 16, cex = 0.55, col = cols)
    } else {
      # DP-colored base-R points (not used for GQ-only runs, but kept for completeness)
      dp_levels <- sort(unique(d$dp_thr[!is.na(d$dp_thr)]))
      if (length(dp_levels) > 1) {
        dp_pal <- grDevices::hcl.colors(length(dp_levels), palette = "viridis", alpha = 1)
        dp_pal2 <- grDevices::adjustcolor(dp_pal, alpha.f = 0.75)
        dcols <- dp_pal2[match(d$dp_thr, dp_levels)]
        points(d$recall, d$precision, pch = 16, cex = 0.55, col = dcols)
      } else {
        points(d$recall, d$precision, pch = 16, cex = 0.55, col = rgb(0, 0, 0, alpha = 0.5))
      }
    }

  }
}

base_out <- file.path(out_dir, "vcfeval_grid_pr_by_sample_base.png")
png(base_out, width = 1650, height = 520, res = 150, bg = "white")
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 1, 0))

gq_levels <- sort(unique(long$gq_thr))
pal <- grDevices::hcl.colors(length(gq_levels), palette = "viridis", alpha = 1)
pal2 <- grDevices::adjustcolor(pal, alpha.f = 0.75)

draw_panel(long, "SNV", c(0.95, 1.00), c(0.95, 1.00), "SNV")
draw_panel(long, "Insertion", c(0.65, 1.00), c(0.75, 1.00), "Insertion")
draw_panel(long, "Deletion", c(0.65, 1.00), c(0.75, 1.00), "Deletion")
mtext("Precision vs recall (vcfeval), one translucent line per sample", outer = TRUE, cex = 0.95)
par(op)
dev.off()
message("wrote ", base_out)
