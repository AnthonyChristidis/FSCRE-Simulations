# -----------------------------------
# Generate CPU Scalability Plots
# -----------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# -----------------------------------

# __________________________
# 1. Load and Compile Data
# __________________________

cat("Loading timing results...\n")

# Path to where your timing files are saved
results_dir <- "results"

# Find all files matching the timing pattern
timing_files <- list.files(results_dir, pattern = "^timing_n=.*\\.rds$", full.names = TRUE)

if(length(timing_files) == 0) {
  stop("No timing files found in the results/ directory.")
}

# Read and combine all files into one data frame
timing_list <- lapply(timing_files, readRDS)
full_timing_df <- do.call(rbind, timing_list)

# _____________________
# 2. Data Aggregation
# _____________________

cat("Aggregating data...\n")

# Calculate the median CPU time for each configuration to be robust to cluster noise
summary_df <- full_timing_df %>%
  group_by(n, p, Method) %>%
  summarize(
    Median_Time = median(Time),
    .groups = 'drop'
  )

# Rename methods for the plot legend
summary_df$Method <- factor(summary_df$Method, 
                            levels = c("FSCRE", "DDC_EN"), 
                            labels = c("FSCRE", "DDC-EN"))

# _________________
# 3. Create Plots
# _________________

cat("Generating plots...\n")

# Define custom theme for publication (Black and White, clean)
pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"), # Slightly larger text for readability
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.width = unit(1.5, "cm"), # Make legend lines longer so dashes are visible
    panel.grid.minor = element_blank()
  )

# --- Panel A: Scaling with p (Fixed n = 100) ---
fixed_n <- 100
df_p <- summary_df %>% filter(n == fixed_n)

plot_p <- ggplot(df_p, aes(x = p, y = Median_Time, group = Method, shape = Method, linetype = Method)) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(50, 100, 500, 1000, 5000)) +
  scale_y_log10() +
  scale_shape_manual(values = c(16, 15)) +     # 16 = Solid Circle, 15 = Solid Square
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    x = expression(bold("Number of Predictors (") * italic(p) * bold(")")), # Math formatting
    y = "Median CPU Time (Seconds)"
  ) +
  pub_theme

# --- Panel B: Scaling with n (Fixed p = 1000) ---
fixed_p <- 1000
df_n <- summary_df %>% filter(p == fixed_p)

plot_n <- ggplot(df_n, aes(x = n, y = Median_Time, group = Method, shape = Method, linetype = Method)) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(50, 100, 200, 500)) +
  scale_y_log10() +
  scale_shape_manual(values = c(16, 15)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    x = expression(bold("Sample Size (") * italic(n) * bold(")")), # Math formatting
    y = "Median CPU Time (Seconds)"
  ) +
  pub_theme

# _____________________
# 4. Combine and Save
# _____________________

cat("Saving combined plot...\n")

# Extract the legend from one plot
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(plot_p + theme(legend.position="bottom"))

# Arrange the plots side-by-side in the top row, and the legend in the bottom row
combined_plot <- grid.arrange(
  arrangeGrob(plot_p + theme(legend.position="none"), 
              plot_n + theme(legend.position="none"), 
              nrow=1),
  mylegend, 
  nrow=2, 
  heights=c(10, 1) 
)

# Save as PDF 
ggsave("figures/CPU_Scalability_Plot.pdf", plot = combined_plot, width = 10, height = 5, units = "in")

cat("Done. Plots saved in the results/ directory.\n")