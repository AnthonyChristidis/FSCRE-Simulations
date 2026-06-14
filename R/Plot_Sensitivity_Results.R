# -------------------------------------------------------------
# Plot Sensitivity Analysis for FSCRE Number of Models (K)
# -------------------------------------------------------------

# Clear all memory
rm(list = ls())

# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# -------------------------------------------------------------

# _____________________________________________________
# 1. Process Results and Generate Presentation Figure
# _____________________________________________________

# Parameters used in the simulation (to find the files)
scenario_val <- "mixture_correlation"
snr_val <- 1.0
contam_str <- "0.1_0.05"
p_active_vec <- c(50, 100, 200)
K_vec <- 1:20

# Initialize an empty data frame to store aggregated results
summary_data <- data.frame()

# Loop through the saved files and extract median performance
for (p_active_val in p_active_vec) {
  
  filename <- paste0("results/sensitivity_K_scen=", scenario_val, 
                     "_snr=", snr_val, 
                     "_pAct=", p_active_val,
                     "_contam=", contam_str, ".rds")
  
  if (file.exists(filename)) {
    res_array <- readRDS(filename)
    res_median <- apply(res_array, c(1, 2), median, na.rm = TRUE)
    
    df_temp <- as.data.frame(res_median)
    df_temp$K <- K_vec
    df_temp$p_active <- paste0("Active Predictors: ", p_active_val, " (", (p_active_val/500)*100, "%)")
    
    summary_data <- bind_rows(summary_data, df_temp)
  } else {
    warning(paste("File not found:", filename))
  }
}

# Check if data was successfully loaded
if (nrow(summary_data) == 0) {
  stop("No results data found. Please run Generate_Sensitivity_Results.R first.")
}

# _____________________________
# 2. Data Formatting for ggplot2
# _____________________________

plot_data <- summary_data |>
  select(K, p_active, MSPE, RC, PR) |>
  pivot_longer(cols = c(MSPE, RC, PR), names_to = "Metric", values_to = "Value")

plot_data$Metric <- factor(plot_data$Metric, 
                           levels = c("MSPE", "RC", "PR"), 
                           labels = c("Prediction Error (MSPE)", "Recall", "Precision"))

plot_data$p_active <- factor(plot_data$p_active, 
                             levels = c("Active Predictors: 50 (10%)", 
                                        "Active Predictors: 100 (20%)", 
                                        "Active Predictors: 200 (40%)"))

# Extract the RLARS Baseline (K=1) to plot as a horizontal line
rlars_data <- plot_data |> filter(K == 1)

# 
# --- PUBLICATION-READY THEME & COLORS ---
# 

pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12, color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  )

fscre_color <- "#0072B2"
rlars_color <- "#D55E00"

# _____________________
# 3. Generate the Figure
# _____________________

if (!dir.exists("figures")) dir.create("figures")

p <- ggplot(plot_data, aes(x = K, y = Value)) +
  # Add the RLARS (K=1) horizontal baseline
  geom_hline(data = rlars_data, aes(yintercept = Value, linetype = "RLARS (Single Model)"), 
             color = rlars_color, linewidth = 1) +
  # Add the FSCRE curve
  geom_line(aes(color = "FSCRE (Ensemble)"), linewidth = 1.2) +
  geom_point(aes(color = "FSCRE (Ensemble)"), size = 2.5) +
  facet_grid(Metric ~ p_active, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_color_manual(values = c("FSCRE (Ensemble)" = fscre_color)) +
  scale_linetype_manual(values = c("RLARS (Single Model)" = "dashed")) +
  pub_theme +
  labs(
    x = expression(bold("Number of Models (") * italic(K) * bold(")")),
    y = "Median Performance Value"
  )

print(p)

ggsave(filename = "figures/FSCRE_K_Sensitivity.pdf", plot = p, 
       width = 10, height = 7, units = "in", dpi = 300)

cat("\nFigure saved to: figures/FSCRE_K_Sensitivity.pdf\n")