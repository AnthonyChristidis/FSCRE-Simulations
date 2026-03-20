# ----------------------------------------------------
# Generate Bioinformatics Application Plots (MSPE Boxplots)
# ----------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------------------

# __________________________
# 1. Load Data
# __________________________

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

# Helper function to melt matrix into tidy dataframe
melt_mspe <- function(mat) {
  df <- as.data.frame(mat)
  df$Split <- 1:nrow(df)
  
  # Pivot longer
  df_long <- pivot_longer(df, 
                          cols = -Split, 
                          names_to = "Condition_Method", 
                          values_to = "MSPE")
  
  # Separate the "Orig/Contam" prefix from the "Method" name
  df_long <- df_long %>%
    separate(Condition_Method, into = c("Condition", "Method"), sep = "_", extra = "merge")
  
  return(df_long)
}

plot_data <- melt_mspe(er_res)

# _________________________________________
# 3. Filter and Format Data for Plotting
# _________________________________________

# Remove the methods with astronomical errors so the plot is readable
# Sparse_S and CR_Lasso are excluded
methods_to_keep <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "RLARS", "FSCRE")
plot_data <- plot_data %>% filter(Method %in% methods_to_keep)

# Clean up names for the plot labels
plot_data$Method <- factor(plot_data$Method, 
                           levels = c("ElasticNet", "DDC_EN", "DDC_RGLM", "RLARS", "FSCRE"),
                           labels = c("EN", "DDC-EN", "DDC-RGLM", "RLARS", "FSCRE"))

plot_data$Condition <- factor(plot_data$Condition, 
                              levels = c("Orig", "Contam"),
                              labels = c("Original Data", "Contaminated Data"))

# _________________
# 4. Generate Plot
# _________________

cat("Generating plot...\n")

# Define a clean, professional black and white theme
pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank() 
  )

final_plot <- ggplot(plot_data, aes(x = Method, y = MSPE, fill = Condition)) +
  geom_boxplot(width = 0.6, outlier.size = 1.5, outlier.alpha = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Original Data" = "grey80", "Contaminated Data" = "grey40")) +
  labs(
    x = NULL,
    y = "MSPE"
  ) +
  pub_theme

# __________________
# 5. Save the Plot
# __________________

cat("Saving plot...\n")

# Ensure figures directory exists
if (!dir.exists("figures")) dir.create("figures")

# Save plot
ggsave("figures/Application_MSPE_Plot.pdf", plot = final_plot, width = 7, height = 6, units = "in")

cat("Done. Plots saved to figures/ directory.\n")

# ______________________________________
# 6. Gene Selection Stability Analysis
# ______________________________________

cat("\n\n--- Gene Selection Stability Analysis (ER-alpha) ---\n")

# Path to the genes file (assuming you ran the Application2 script)
genes_file <- "results/Application_TCGA_ER_alpha_Genes.rds"

if(!file.exists(genes_file)) {
  stop("Missing Genes result file in results/ directory. Did the simulation finish saving it?")
}

genes_res <- readRDS(genes_file)
# genes_res is a list of length 50 (one per split).
# Each element has $Orig and $Contam.
# Each of those has $ElasticNet, $DDC_EN, $DDC_RGLM, $RLARS, $FSCRE (or whichever methods you saved).

# Number of splits
n_splits <- length(genes_res)

# Helper function to count selection frequencies for a specific method and condition
get_selection_props <- function(results_list, condition, method) {
  # Extract the character vectors of selected genes across all 50 splits
  all_genes_selected <- lapply(results_list, function(split) {
    # Check if the split has data for this method (handles potential NAs/crashes in one split)
    if (!is.null(split[[condition]][[method]])) {
      return(split[[condition]][[method]])
    } else {
      return(character(0))
    }
  })
  
  # Flatten the list into a single long vector of all selections
  flat_genes <- unlist(all_genes_selected)
  
  # Count frequencies and convert to proportions
  gene_counts <- table(flat_genes)
  gene_props <- gene_counts / n_splits
  
  return(gene_props)
}

# 1. Get proportions for FSCRE on Contaminated Data (Our primary target)
props_fscre_contam <- get_selection_props(genes_res, "Contam", "FSCRE")

# 2. Get proportions for Elastic Net on Contaminated Data (The baseline)
props_en_contam <- get_selection_props(genes_res, "Contam", "ElasticNet")

# 3. Get proportions for Elastic Net on Original Data (The "Clean" baseline)
props_en_orig <- get_selection_props(genes_res, "Orig", "ElasticNet")

# 4. Get proportions for DDC+EN on Contaminated Data (The "Simple Robust" baseline)
props_ddc_en_contam <- get_selection_props(genes_res, "Contam", "DDC_EN")


# --- Combine into a single comparison table ---

all_fscre_genes <- names(props_fscre_contam)

selection_table <- data.frame(
  Gene = all_fscre_genes,
  FSCRE_Contam = as.numeric(props_fscre_contam),
  stringsAsFactors = FALSE
)

# Add DDC-EN Contam
selection_table$DDC_EN_Contam <- as.numeric(props_ddc_en_contam[selection_table$Gene])
selection_table$DDC_EN_Contam[is.na(selection_table$DDC_EN_Contam)] <- 0

# Add EN Contam
selection_table$EN_Contam <- as.numeric(props_en_contam[selection_table$Gene])
selection_table$EN_Contam[is.na(selection_table$EN_Contam)] <- 0

# Add EN Orig
selection_table$EN_Orig <- as.numeric(props_en_orig[selection_table$Gene])
selection_table$EN_Orig[is.na(selection_table$EN_Orig)] <- 0

# Sort the table descending based on FSCRE
selection_table <- selection_table[order(selection_table$FSCRE_Contam, decreasing = TRUE), ]

cat("\nTop 20 Genes Selected by FSCRE under Targeted Contamination:\n")
print(head(selection_table, 20), row.names = FALSE)

