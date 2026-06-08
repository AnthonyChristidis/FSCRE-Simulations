# ------------------------------------------------------------
# Generate Bioinformatics Application Plots (MSPE Boxplots)
# ------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------

# ______________
# 1. Load Data
# ______________

cat("\n--- Loading MSPE results ---\n")

er_file <- "results/Application_TCGA_ER_alpha_MSPE.rds"

if(!file.exists(er_file)) {
  stop("Missing ER-alpha MSPE result file in results/ directory.")
}

er_res <- readRDS(er_file)

# __________________________
# 2. Reshape Data for Plot
# __________________________

cat("Reshaping data...\n")

melt_mspe <- function(mat) {
  df <- as.data.frame(mat)
  df$Split <- 1:nrow(df)
  
  df_long <- pivot_longer(df, 
                          cols = -Split, 
                          names_to = "Condition_Method", 
                          values_to = "MSPE")
  
  df_long <- df_long %>%
    separate(Condition_Method, into = c("Condition", "Method"), sep = "_", extra = "merge")
  
  return(df_long)
}

plot_data <- melt_mspe(er_res)

# _________________________________________
# 3. Filter and Format Data for Plotting
# _________________________________________

methods_to_keep <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "RLARS", "FSCRE")
plot_data <- plot_data %>% filter(Method %in% methods_to_keep)

plot_data$Method <- factor(plot_data$Method, 
                           levels = c("ElasticNet", "DDC_EN", "DDC_RGLM", "RLARS", "FSCRE"),
                           labels = c("EN", "DDC-EN", "DDC-RGLM", "RLARS", "FSCRE"))

plot_data$Condition <- factor(plot_data$Condition, 
                              levels = c("Orig", "Contam"),
                              labels = c("Original Data", "Targeted Contamination"))

# --- Themes and colors ---

pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13, face = "bold"),
    legend.key.size = unit(1.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed") 
  )

# The "Baseline vs. Stress" Palette
# Grey = Clean baseline. Crimson Red = The stress test (zero orange undertones). 
fill_colors <- c(
  "Original Data"          = "#E0E0E0",  # Light, clean grey
  "Targeted Contamination" = "#C41E3A"   # True Crimson / Cardinal Red
)

# _________________
# 4. Generate Plot
# _________________

cat("Generating plot...\n")

final_plot <- ggplot(plot_data, aes(x = Method, y = MSPE, fill = Condition)) +
  # Setting both widths to 0.5 makes them skinny but perfectly touching
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.5), 
               outlier.size = 1.5, outlier.shape = 21, outlier.alpha = 0.6, 
               alpha = 0.9, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = fill_colors) +
  labs(
    x = NULL,
    y = "Mean Squared Prediction Error (MSPE)"
  ) +
  pub_theme

# __________________
# 5. Save the Plot
# __________________

cat("Saving plot...\n")

if (!dir.exists("figures")) dir.create("figures")

ggsave("figures/Application_MSPE_Plot.pdf", plot = final_plot, width = 8, height = 6, units = "in")

cat("Done. Plots saved to figures/ directory.\n")

# ______________________________________
# 6. Gene Selection Stability Analysis
# ______________________________________

cat("\n\n--- Gene Selection Stability Analysis (ER-alpha) ---\n")

genes_file <- "results/Application_TCGA_ER_alpha_Genes.rds"

if(!file.exists(genes_file)) {
  stop("Missing Genes result file in results/ directory.")
}

genes_res <- readRDS(genes_file)
n_splits <- length(genes_res)

get_selection_props <- function(results_list, condition, method) {
  all_genes_selected <- lapply(results_list, function(split) {
    if (!is.null(split[[condition]][[method]])) {
      return(split[[condition]][[method]])
    } else {
      return(character(0))
    }
  })
  
  flat_genes <- unlist(all_genes_selected)
  gene_counts <- table(flat_genes)
  gene_props <- gene_counts / n_splits
  return(gene_props)
}

props_fscre_contam <- get_selection_props(genes_res, "Contam", "FSCRE")
props_en_contam <- get_selection_props(genes_res, "Contam", "ElasticNet")
props_en_orig <- get_selection_props(genes_res, "Orig", "ElasticNet")
props_ddc_en_contam <- get_selection_props(genes_res, "Contam", "DDC_EN")

all_fscre_genes <- names(props_fscre_contam)

selection_table <- data.frame(
  Gene = all_fscre_genes,
  FSCRE_Contam = as.numeric(props_fscre_contam),
  stringsAsFactors = FALSE
)

selection_table$DDC_EN_Contam <- as.numeric(props_ddc_en_contam[selection_table$Gene])
selection_table$DDC_EN_Contam[is.na(selection_table$DDC_EN_Contam)] <- 0

selection_table$EN_Contam <- as.numeric(props_en_contam[selection_table$Gene])
selection_table$EN_Contam[is.na(selection_table$EN_Contam)] <- 0

selection_table$EN_Orig <- as.numeric(props_en_orig[selection_table$Gene])
selection_table$EN_Orig[is.na(selection_table$EN_Orig)] <- 0

selection_table <- selection_table[order(selection_table$FSCRE_Contam, decreasing = TRUE), ]

cat("\nTop 20 Genes Selected by FSCRE under Targeted Contamination:\n")
print(head(selection_table, 20), row.names = FALSE)