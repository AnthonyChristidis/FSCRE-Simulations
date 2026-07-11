# --------------------------------------------------------------
# Generate Simulation Plots (MSPE, Precision, Recall Grids)
# --------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# --------------------------------------------------------------

# __________________________
# 1. Load and Compile Data
# __________________________

cat("\n--- Loading simulation results across SNRs and Sparsity ---\n")

snrs <- c(0.5, 1, 2)
p_actives <- c(50, 100, 200)
scenario <- "mixture_correlation"
contam_str <- "0.1_0.05"

# Load ALL relevant methods so we can filter them specifically for different plots
# Keeping "RLARS" here to extract data correctly from the RDS files
methods_to_load <- c("DDC_EN", "DDC_RGLM", "RLARS", "FSCRE")
df_list <- list()

for (s in snrs) {
  for (pa in p_actives) {
    
    file_name <- sprintf("results/res_scen=%s_snr=%s_pAct=%s_contam=%s.rds", 
                         scenario, s, pa, contam_str)
    
    if(!file.exists(file_name)) {
      warning(sprintf("Missing file: %s. Skipping...", file_name))
      next
    }
    
    res_array <- readRDS(file_name)
    reps <- dim(res_array)[3]
    
    for (m in methods_to_load) {
      if (m %in% rownames(res_array)) {
        temp_df <- data.frame(
          Method = rep(m, reps),
          Rep = 1:reps,
          SNR = s,
          Sparsity = pa,
          MSPE = res_array[m, "MSPE", ],
          Precision = res_array[m, "PR", ],
          Recall = res_array[m, "RC", ]
        )
        df_list[[length(df_list) + 1]] <- temp_df
      }
    }
  }
}

plot_data <- do.call(rbind, df_list)

# _____________________________
# 2. Format Data for Plotting
# _____________________________

cat("Formatting data...\n")

# Standardize method names. Notice "RLARS" is mapped to the label "CellRLARS" here.
plot_data$Method <- factor(plot_data$Method, 
                           levels = c("DDC_EN", "DDC_RGLM", "RLARS", "FSCRE"),
                           labels = c("DDC-EN", "DDC-RGLM", "CellRLARS", "FSCRE"))

# Create nice labels for the facets
plot_data$SNR_Label <- factor(plot_data$SNR,
                              levels = c(0.5, 1, 2),
                              labels = c("Low Signal (SNR = 0.5)", 
                                         "Moderate Signal (SNR = 1.0)", 
                                         "High Signal (SNR = 2.0)"))

# Ensure sparsity is treated as a discrete category for the X-axis
plot_data$Sparsity_Label <- factor(plot_data$Sparsity, levels = c(50, 100, 200))

# --- Themes and colors ---

# Define a clean, highly professional black and white theme
pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 13), 
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed") 
  )

# Professional, colorblind-friendly palette (Okabe-Ito inspired)
# Changed "RLARS" key to "CellRLARS" to match the factor labels
fill_colors <- c(
  "DDC-EN"    = "#D55E00",  # Rust / Vermilion (Strong baseline contrast)
  "DDC-RGLM"  = "#E69F00",  # Golden Yellow (Secondary baseline)
  "CellRLARS" = "#56B4E9",  # Sky Blue (Links it visually to the proposed method)
  "FSCRE"     = "#0072B2"   # Deep Yale Blue (Authoritative, stable proposed method)
)

# _________________
# 3. Generate Plots
# _________________

cat("Generating plots...\n")

if (!dir.exists("figures")) dir.create("figures")

# --- FIGURE 1: MSPE Plot (1x3 Facets) ---

# Filter specifically for the MSPE methods (excluding CellRLARS to avoid clutter)
df_mspe <- plot_data |> filter(Method %in% c("DDC-EN", "DDC-RGLM", "FSCRE"))

# Calculate dynamic y-axis limits to trim ONLY the most extreme, scale-breaking outliers.
max_y_limit <- quantile(df_mspe$MSPE, 0.99, na.rm = TRUE) * 1.1
min_y_limit <- min(df_mspe$MSPE, na.rm = TRUE) * 0.90

p_mspe <- ggplot(df_mspe, aes(x = Sparsity_Label, y = MSPE, fill = Method)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, outlier.alpha = 0.6, outlier.shape = 21, alpha = 0.9, color = "black", linewidth = 0.6) +
  facet_grid(. ~ SNR_Label, scales = "free_y") + 
  scale_fill_manual(values = fill_colors) +
  coord_cartesian(ylim = c(min_y_limit, max_y_limit)) + 
  labs(
    x = "Number of Active Predictors",
    y = "Mean Squared Prediction Error (MSPE)"
  ) +
  pub_theme

ggsave("figures/Simulation_MSPE_Grid.pdf", plot = p_mspe, width = 12, height = 5, units = "in")

# --- FIGURE 2: Recall and Precision Plot (2x3 Facets) ---

# Filter to include DDC-EN, CellRLARS, and FSCRE
df_rcpr <- plot_data |> filter(Method %in% c("DDC-EN", "CellRLARS", "FSCRE"))

# Melt the data so RC and PR are in the same column for faceting
df_rcpr_long <- df_rcpr |>
  pivot_longer(cols = c(Recall, Precision), names_to = "Metric", values_to = "Value")

# Order the metric factor so Recall is top row, Precision is bottom row
df_rcpr_long$Metric <- factor(df_rcpr_long$Metric, levels = c("Recall", "Precision"))

p_rcpr <- ggplot(df_rcpr_long, aes(x = Sparsity_Label, y = Value, fill = Method)) +
  geom_boxplot(width = 0.6, outlier.size = 1.5, outlier.alpha = 0.6, outlier.shape = 21, alpha = 0.9, color = "black", linewidth = 0.6) +
  facet_grid(Metric ~ SNR_Label, scales = "free_y") + 
  scale_fill_manual(values = fill_colors) +
  labs(
    x = "Number of Active Predictors",
    y = NULL
  ) +
  pub_theme

ggsave("figures/Simulation_RCPR_Grid.pdf", plot = p_rcpr, width = 12, height = 7, units = "in")

cat("Done. Plots saved to figures/ directory.\n")