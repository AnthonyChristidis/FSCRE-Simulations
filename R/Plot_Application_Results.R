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
  
  df_long <- df_long |>
    separate(Condition_Method, into = c("Condition", "Method"), sep = "_", extra = "merge")
  
  return(df_long)
}

plot_data <- melt_mspe(er_res)

# _________________________________________
# 3. Filter and Format Data for Plotting
# _________________________________________

# We keep "RLARS" here so it can extract from the RDS properly
methods_to_keep <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "RLARS", "FSCRE")
plot_data <- plot_data |> filter(Method %in% methods_to_keep)

# We update the 'labels' argument to change "RLARS" to "CellRLARS" in the plot
plot_data$Method <- factor(plot_data$Method, 
                           levels = c("ElasticNet", "DDC_EN", "DDC_RGLM", "RLARS", "FSCRE"),
                           labels = c("EN", "DDC-EN", "DDC-RGLM", "CellRLARS", "FSCRE"))

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

# Methods you want to summarize (keeping "RLARS" internally to match the RDS file keys)
methods_all <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "Sparse_S", "CR_Lasso", "RLARS", "FSCRE")
conditions_all <- c("Orig", "Contam")

get_selection_props <- function(results_list, condition, method) {
  all_genes_selected <- lapply(results_list, function(split) {
    # split is like: list(Orig = res_orig$vars, Contam = res_cont$vars)
    if (!is.null(split[[condition]]) &&
        !is.null(split[[condition]][[method]]) &&
        length(split[[condition]][[method]]) > 0) {
      return(split[[condition]][[method]])
    } else {
      return(character(0))
    }
  })

  flat_genes <- unlist(all_genes_selected, use.names = FALSE)

  if (length(flat_genes) == 0) {
    return(setNames(numeric(0), character(0)))
  }

  gene_counts <- table(flat_genes)
  gene_props <- gene_counts / n_splits
  return(gene_props)
}

# 1) Compute proportions for every (condition, method)
props_list <- list()
for (cond in conditions_all) {
  for (meth in methods_all) {
    nm <- paste0(meth, "_", cond)
    props_list[[nm]] <- get_selection_props(genes_res, cond, meth)
  }
}

# 2) Union of all genes that were ever selected by any method in any condition
all_genes <- sort(unique(unlist(lapply(props_list, names), use.names = FALSE)))

# 3) Build wide table: rows=genes, cols=Method_Orig / Method_Contam
selection_table <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)

for (nm in names(props_list)) {
  v <- props_list[[nm]]
  selection_table[[nm]] <- as.numeric(v[selection_table$Gene])
  selection_table[[nm]][is.na(selection_table[[nm]])] <- 0
}

# Rename RLARS columns to CellRLARS to match the manuscript
names(selection_table) <- gsub("RLARS", "CellRLARS", names(selection_table))

# Optional: add a simple summary column (max selection rate across all columns)
selection_table$MaxProp_Any <- apply(selection_table[, setdiff(names(selection_table), "Gene")], 1, max)

# Optional: filter to keep table manageable (tweak threshold if desired)
# selection_table <- subset(selection_table, MaxProp_Any >= 0.10)

# 4) Sort & print: choose which column drives the ranking
rank_col <- "FSCRE_Contam"   # change this if you want, e.g. "ElasticNet_Orig"
if (!(rank_col %in% names(selection_table))) {
  stop("rank_col not found in selection_table: ", rank_col)
}

selection_table <- selection_table[order(selection_table[[rank_col]], decreasing = TRUE), ]

cat("\nTop 20 genes by ", rank_col, " (with Orig/Contam proportions for all methods):\n", sep = "")
head(selection_table[, c("Gene",
                          "FSCRE_Contam","FSCRE_Orig",
                          "CellRLARS_Contam","CellRLARS_Orig",
                          "ElasticNet_Contam","ElasticNet_Orig",
                          "DDC_EN_Contam","DDC_EN_Orig",
                          "DDC_RGLM_Contam","DDC_RGLM_Orig",
                          "Sparse_S_Contam","Sparse_S_Orig",
                          "CR_Lasso_Contam","CR_Lasso_Orig")],
      row.names = FALSE, n = 20)