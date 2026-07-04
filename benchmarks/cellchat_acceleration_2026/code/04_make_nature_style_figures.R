#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(grid)
})

root <- Sys.getenv("ROOT", "/home/dzf/cellchat_acceleration")
tab_dir <- file.path(root, "results", "tables")
fig_dir <- file.path(root, "results", "figures_nature")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

m <- fread(file.path(tab_dir, "all_metrics.tsv"), na.strings = c("NA", "", "NaN"))
scale_levels <- c("1000", "5000", "10000", "25000", "50000", "all")
scale_labels <- c("1k", "5k", "10k", "25k", "50k", "all")
m[, scale := factor(as.character(n_cells_target), levels = scale_levels, labels = scale_labels)]
m[, scale_rev := factor(as.character(n_cells_target), levels = rev(scale_levels), labels = rev(scale_labels))]

ok <- m[status %in% c("ok", "ok_resumed")]
full <- ok[ablation == "full" & engine %in% c("baseline", "accelerated")]
wide <- dcast(full, experiment_id + dataset_id + n_cells_target + scale + scale_rev + seed + metrics_file ~ engine, value.var = "elapsed_sec")
wide <- wide[!is.na(baseline) & !is.na(accelerated)]
wide[, speedup := baseline / accelerated]

scale_sum <- wide[, .(
  baseline_med = median(baseline),
  baseline_q25 = quantile(baseline, 0.25),
  baseline_q75 = quantile(baseline, 0.75),
  accel_med = median(accelerated),
  accel_q25 = quantile(accelerated, 0.25),
  accel_q75 = quantile(accelerated, 0.75),
  speedup_med = median(speedup),
  speedup_q25 = quantile(speedup, 0.25),
  speedup_q75 = quantile(speedup, 0.75),
  n = .N
), by = .(n_cells_target, scale, scale_rev)]
setorder(scale_sum, scale)

runtime_long <- rbindlist(list(
  scale_sum[, .(scale, scale_rev, engine = "Original CellChat", median = baseline_med, q25 = baseline_q25, q75 = baseline_q75)],
  scale_sum[, .(scale, scale_rev, engine = "CellChatAccelRcpp", median = accel_med, q25 = accel_q25, q75 = accel_q75)]
))
runtime_long[, engine := factor(engine, levels = c("Original CellChat", "CellChatAccelRcpp"))]

acc <- ok[engine == "comparison" & ablation == "full"]
full_accel <- full[engine == "accelerated", .(dataset_id, n_cells_target, scale, seed, full_elapsed = elapsed_sec)]
abl <- ok[engine == "accelerated" & ablation %in% c("no_accel_kernel", "no_accel_pathway", "no_accel_aggregate"),
          .(dataset_id, n_cells_target, scale, seed, ablation, elapsed_sec)]
abl <- merge(abl, full_accel, by = c("dataset_id", "n_cells_target", "scale", "seed"), all.x = TRUE)
abl[, slowdown := elapsed_sec / full_elapsed]
abl[, component := factor(ablation,
  levels = c("no_accel_kernel", "no_accel_pathway", "no_accel_aggregate"),
  labels = c("Probability kernel", "Pathway aggregation", "Network aggregation")
)]
abl_summary <- abl[!is.na(slowdown), .(
  median = median(slowdown),
  q25 = quantile(slowdown, 0.25),
  q75 = quantile(slowdown, 0.75),
  n = .N
), by = component]

heat <- wide[, .(median_speedup = median(speedup, na.rm = TRUE)), by = .(dataset_id, scale)]
heat[, dataset_short := sub("^3CA_data__", "", dataset_id)]
heat[, dataset_short := sub("^normal_control__", "", dataset_short)]
ord <- heat[, .(score = median(median_speedup, na.rm = TRUE)), by = dataset_short][order(score)]$dataset_short
heat[, dataset_short := factor(dataset_short, levels = ord)]
heat[, speed_label := sprintf("%.0f", median_speedup)]

paired_n <- nrow(wide)
overall_speedup <- median(wide$speedup)
max_diff <- max(acc$max_abs_prob_diff, na.rm = TRUE)
min_pearson <- min(acc$pearson_prob, na.rm = TRUE)
abl_ymax <- max(abl$slowdown, na.rm = TRUE) * 1.18

# Muted red-blue palette for publication figures: original CellChat in red,
# accelerated CellChatAccelRcpp in blue; ordinal scales use warm-to-cool tones.
pal_engine <- c("Original CellChat" = "#B54A4A", "CellChatAccelRcpp" = "#2F6DB3")
pal_scale <- c("1k" = "#B54A4A", "5k" = "#D27467", "10k" = "#E9B0A5", "25k" = "#AFCBE8", "50k" = "#6E9ED0", "all" = "#2F5F9E")
pal_abl <- c("Probability kernel" = "#B54A4A", "Pathway aggregation" = "#6E7EBE", "Network aggregation" = "#2F6DB3")

nature_theme <- theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 11.2, colour = "#111111", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 8.5, colour = "#555555", margin = margin(b = 6)),
    axis.title = element_text(size = 9.2, colour = "#222222"),
    axis.text = element_text(size = 8.2, colour = "#333333"),
    axis.line = element_line(linewidth = 0.35, colour = "#303030"),
    axis.ticks = element_line(linewidth = 0.28, colour = "#303030"),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "#ECECEC"),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.2),
    legend.key.width = unit(8, "mm"),
    legend.key.height = unit(3.6, "mm"),
    plot.margin = margin(7, 9, 7, 9)
  )

save_both <- function(plot, name, width, height) {
  pdf <- file.path(fig_dir, paste0(name, ".pdf"))
  png <- file.path(fig_dir, paste0(name, ".png"))
  ggsave(pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png, plot, width = width, height = height, units = "in", dpi = 500, bg = "white")
  cat(pdf, "\n", png, "\n", sep = "")
}

# Figure 1: runtime compression dumbbell
fig1 <- ggplot(scale_sum, aes(y = scale_rev)) +
  geom_segment(aes(x = accel_med, xend = baseline_med, yend = scale_rev), colour = "#D7D3D0", linewidth = 2.6, lineend = "round") +
  geom_errorbarh(data = runtime_long, aes(xmin = q25, xmax = q75, y = scale_rev, colour = engine), height = 0.13, linewidth = 0.52) +
  geom_point(data = runtime_long, aes(x = median, y = scale_rev, colour = engine), size = 3.1) +
  geom_text(aes(x = 1850, label = sprintf("%.1fx", speedup_med)), hjust = 1, size = 3.05, fontface = "bold", colour = "#222222") +
  annotate("text", x = 1850, y = 6.25, label = "median speedup", hjust = 1, size = 2.75, colour = "#6A6A6A") +
  scale_x_log10(labels = label_number(), breaks = c(3, 10, 30, 100, 300, 1000, 2000), limits = c(2.5, 2200)) +
  scale_colour_manual(values = pal_engine) +
  labs(
    title = "Runtime compression across cell scales",
    subtitle = sprintf("%d paired original/accelerated CellChat runs; medians with interquartile ranges.", paired_n),
    x = "Elapsed time (s, log10)", y = "Target cells"
  ) +
  nature_theme

# Figure 2: speedup distribution
fig2 <- ggplot(wide, aes(scale, speedup, fill = scale)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35, colour = "#9A9A9A") +
  geom_violin(width = 0.82, colour = NA, alpha = 0.93, trim = FALSE) +
  geom_boxplot(width = 0.17, outlier.shape = NA, fill = "white", linewidth = 0.35, colour = "#2B2F38") +
  geom_point(position = position_jitter(width = 0.055, height = 0), size = 0.85, alpha = 0.36, colour = "#202020") +
  geom_text(data = scale_sum, aes(x = scale, y = speedup_q75 * 1.45, label = sprintf("%.1fx", speedup_med)), inherit.aes = FALSE, size = 3.0, fontface = "bold", colour = "#1A1A1A") +
  scale_y_log10(labels = label_number(accuracy = 1), breaks = c(1, 3, 10, 30, 60), limits = c(0.9, 75)) +
  scale_fill_manual(values = pal_scale) +
  labs(
    title = "Speedup distribution",
    subtitle = sprintf("Overall median speedup %.1fx; boxes show median and interquartile range.", overall_speedup),
    x = "Target cells", y = "Original / accelerated runtime"
  ) +
  nature_theme + theme(legend.position = "none", panel.grid.major.x = element_blank())

# Figure 3: paired runtime scatter
contours <- data.table(mult = c(1, 5, 10, 30), label = c("1x", "5x", "10x", "30x"))
fig3 <- ggplot(wide, aes(baseline, accelerated, colour = scale)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35, colour = "#888888") +
  geom_abline(data = contours[mult > 1], aes(slope = 1 / mult, intercept = 0), inherit.aes = FALSE, linetype = 3, linewidth = 0.28, colour = "#B7B7B7") +
  annotate("text", x = c(35, 120, 360), y = c(7, 12, 12), label = c("5x", "10x", "30x"), size = 2.8, colour = "#777777") +
  geom_point(size = 2.05, alpha = 0.78) +
  scale_x_log10(labels = label_number(), breaks = c(30, 100, 300, 1000, 3000)) +
  scale_y_log10(labels = label_number(), breaks = c(3, 10, 30, 100, 300)) +
  scale_colour_manual(values = pal_scale) +
  coord_equal(xlim = c(25, 3500), ylim = c(2.5, 450), expand = TRUE) +
  labs(
    title = "Paired runtime comparison",
    subtitle = "Log-log paired runtimes; dashed guides mark equal runtime and constant speedup.",
    x = "Original CellChat elapsed time (s)", y = "CellChatAccelRcpp elapsed time (s)"
  ) +
  nature_theme + theme(panel.grid.major = element_line(linewidth = 0.22, colour = "#EDEDED"))

# Figure 4: numerical equivalence
fig4 <- ggplot(acc, aes(scale, max_abs_prob_diff, fill = scale)) +
  geom_violin(width = 0.78, colour = NA, alpha = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.17, outlier.shape = NA, fill = "white", linewidth = 0.34, colour = "#2B2F38") +
  geom_point(position = position_jitter(width = 0.05, height = 0), size = 0.72, alpha = 0.38, colour = "#202020") +
  annotate("label", x = 3.72, y = 1.82e-16, label = sprintf("Pearson = %.3f in all paired runs\nmax |delta prob| = %.2e\nall < 1e-12 tolerance", min_pearson, max_diff), label.size = 0.22, fill = "white", colour = "#2B2F38", size = 2.75, hjust = 0) +
  scale_y_log10(labels = label_scientific(), breaks = c(4e-17, 6e-17, 8e-17, 1e-16, 1.4e-16)) +
  coord_cartesian(ylim = c(3.8e-17, 1.65e-16)) +
  scale_fill_manual(values = pal_scale) +
  labs(
    title = "Numerical equivalence to original CellChat",
    subtitle = "Zoomed view of floating-point-scale probability differences across all paired benchmarks.",
    x = "Target cells", y = "Max absolute probability difference"
  ) +
  nature_theme + theme(legend.position = "none")

# Figure 5: component ablation
fig5 <- ggplot(abl[!is.na(slowdown)], aes(component, slowdown, fill = component)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.36, colour = "#9A9A9A") +
  geom_violin(width = 0.82, colour = NA, alpha = 0.90, trim = FALSE) +
  geom_boxplot(width = 0.17, outlier.shape = NA, fill = "white", linewidth = 0.35, colour = "#2B2F38") +
  geom_point(aes(colour = scale), position = position_jitter(width = 0.055, height = 0), size = 0.68, alpha = 0.28) +
  geom_text(data = abl_summary, aes(component, pmin(q75 * 1.35, abl_ymax / 1.25), label = sprintf("%.1fx", median)), inherit.aes = FALSE, size = 3.0, fontface = "bold", colour = "#222222") +
  scale_y_log10(labels = label_number(accuracy = 0.1), breaks = c(1, 3, 10, 30, 100)) +
  coord_cartesian(ylim = c(0.85, abl_ymax)) +
  scale_fill_manual(values = pal_abl) +
  scale_colour_manual(values = pal_scale) +
  labs(
    title = "Component ablation identifies the dominant acceleration target",
    subtitle = "Runtime after disabling each accelerated component, normalized to the full accelerated path.",
    x = NULL, y = "Runtime relative to full accelerated path"
  ) +
  nature_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 12, hjust = 1))

# Figure 6: dataset speedup heatmap
fig6 <- ggplot(heat, aes(scale, dataset_short, fill = median_speedup)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = speed_label), size = 2.75, colour = "#1F252C", fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#F8E6E1", "#E7A79D", "#D16C62", "#B8CDE8", "#6E9ED0", "#2F5F9E"),
    trans = "log10", labels = label_number(accuracy = 1), breaks = c(5, 10, 20, 40, 60)
  ) +
  labs(
    title = "Dataset-level speedup landscape",
    subtitle = "Tile values are median speedups across three repeats.",
    x = "Target cells", y = NULL, fill = "Speedup"
  ) +
  nature_theme + theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 7.8),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

save_both(fig1, "Fig01_runtime_compression", 6.6, 4.5)
save_both(fig2, "Fig02_speedup_distribution", 5.8, 4.5)
save_both(fig3, "Fig03_paired_runtime_scatter", 5.7, 4.8)
save_both(fig4, "Fig04_numerical_equivalence", 5.8, 4.5)
save_both(fig5, "Fig05_component_ablation", 6.0, 4.5)
save_both(fig6, "Fig06_dataset_speedup_heatmap", 6.2, 5.5)

cat("figures_dir=", fig_dir, "\n", sep = "")
cat("paired_comparisons=", paired_n, "\n", sep = "")
cat("median_speedup=", sprintf("%.2f", overall_speedup), "\n", sep = "")
cat("max_abs_prob_diff=", sprintf("%.3e", max_diff), "\n", sep = "")
