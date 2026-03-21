# ----------------------------------------------
# Generate Simulation Tables (Average Ranks)
# ----------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)

# ----------------------------------------------

# ___________________________________________________________
# 1. Helper: Compute Ranks for a Single Configuration Array
# ___________________________________________________________

computeRanksForArray <- function(res_array) {
  # res_array dimensions: [Methods, Metrics(MSPE, RC, PR, CPU), Reps]
  
  # Calculate mean across the 50 replications
  mean_metrics <- apply(res_array, c(1, 2), mean, na.rm = TRUE)
  
  # Initialize rank matrix
  rank_mat <- matrix(NA, nrow = nrow(mean_metrics), ncol = 3)
  rownames(rank_mat) <- rownames(mean_metrics)
  colnames(rank_mat) <- c("MSPE", "RC", "PR")
  
  # MSPE (Lower is better) -> rank(x)
  # All methods compete on MSPE
  if ("MSPE" %in% colnames(mean_metrics)) {
    rank_mat[, "MSPE"] <- rank(mean_metrics[, "MSPE"], na.last = "keep", ties.method = "average")
  }
  
  # RC and PR (Higher is better) -> rank(-x)
  # EXCLUDE DDC_RGLM from ranking for selection metrics
  methods_for_selection <- setdiff(rownames(mean_metrics), "DDC_RGLM")
  
  if ("RC" %in% colnames(mean_metrics)) {
    rc_vals <- mean_metrics[methods_for_selection, "RC"]
    ranks <- rank(-rc_vals, na.last = "keep", ties.method = "average")
    rank_mat[methods_for_selection, "RC"] <- ranks
  }
  
  if ("PR" %in% colnames(mean_metrics)) {
    pr_vals <- mean_metrics[methods_for_selection, "PR"]
    ranks <- rank(-pr_vals, na.last = "keep", ties.method = "average")
    rank_mat[methods_for_selection, "PR"] <- ranks
  }
  
  return(rank_mat)
}

# _____________________________
# 2. Load and Process All Files
# _____________________________

cat("\n--- Processing simulation files and computing ranks ---\n")

results_dir <- "results"
files <- list.files(results_dir, pattern = "^res_scen=.*\\.rds$", full.names = TRUE)

if(length(files) == 0) {
  stop("No simulation result files found.")
}

# List to store rank matrices
all_ranks <- list()

for (f in files) {
  # Parse filename to extract scenario using robust regex
  fname <- basename(f)
  
  scenario_name <- str_match(fname, "scen=(.*?)_snr")[,2]
  contam_val <- str_match(fname, "contam=(.*?)\\.rds")[,2]
  
  # Differentiate Clean vs Casewise
  if (scenario_name == "casewise" && contam_val == "0") {
    broad_scenario <- "Clean"
  } else {
    broad_scenario <- scenario_name
  }
  
  # Load data and compute ranks
  res_array <- readRDS(f)
  rank_mat <- computeRanksForArray(res_array)
  
  # Melt into long format for easy aggregation later
  long_ranks <- as.data.frame(as.table(rank_mat))
  colnames(long_ranks) <- c("Method", "Metric", "Rank")
  long_ranks$Scenario <- broad_scenario
  
  all_ranks[[length(all_ranks) + 1]] <- long_ranks
}

master_rank_df <- do.call(rbind, all_ranks)

# ___________________
# 3. Aggregate Ranks
# ___________________

# Average the ranks over all configurations within each broad scenario
aggregated_ranks <- master_rank_df %>%
  group_by(Scenario, Method, Metric) %>%
  summarize(AvgRank = mean(Rank, na.rm = TRUE), .groups = 'drop')

# _________________________________________
# 4. Generate TABLE 1: Cellwise Scenarios
# _________________________________________

cat("\n\n==========================================================================\n")
cat("TABLE 1: CELLWISE CONTAMINATION (Average Ranks)\n")
cat("==========================================================================\n")

cellwise_scenarios <- c("cellwise_marginal", "cellwise_correlation", "mixture_marginal", "mixture_correlation")
cellwise_methods <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "Sparse_S", "CR_Lasso", "RLARS", "FSCRE")

tab1_data <- aggregated_ranks %>%
  filter(Scenario %in% cellwise_scenarios, Method %in% cellwise_methods) %>%
  pivot_wider(names_from = c(Scenario, Metric), values_from = AvgRank, names_sep = ".") %>%
  # Reorder rows to a logical order
  slice(match(cellwise_methods, Method))

# Reorder columns to group by scenario: MSPE, RC, PR
col_order <- c("Method", 
               paste0("cellwise_marginal.", c("MSPE", "RC", "PR")),
               paste0("cellwise_correlation.", c("MSPE", "RC", "PR")),
               paste0("mixture_marginal.", c("MSPE", "RC", "PR")),
               paste0("mixture_correlation.", c("MSPE", "RC", "PR")))

col_order <- intersect(col_order, colnames(tab1_data))
tab1_data <- tab1_data[, col_order]

# Round to 1 decimal place for printing
tab1_print <- tab1_data %>% mutate(across(where(is.numeric), ~round(., 1)))
print(as.data.frame(tab1_print))

# __________________________________________________
# 5. Generate TABLE 2: Clean and Casewise Scenarios
# __________________________________________________

cat("\n\n==========================================================================\n")
cat("TABLE 2: CLEAN AND CASEWISE SCENARIOS (Average Ranks)\n")
cat("==========================================================================\n")

casewise_scenarios <- c("Clean", "casewise")
casewise_methods <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "Sparse_S", "CR_Lasso", "RLARS", "FSCRE", "PENSE", "SparseLTS")

tab2_data <- aggregated_ranks %>%
  filter(Scenario %in% casewise_scenarios, Method %in% casewise_methods) %>%
  pivot_wider(names_from = c(Scenario, Metric), values_from = AvgRank, names_sep = ".") %>%
  slice(match(casewise_methods, Method))

# Reorder columns
col_order2 <- c("Method", 
               paste0("Clean.", c("MSPE", "RC", "PR")),
               paste0("casewise.", c("MSPE", "RC", "PR")))

col_order2 <- intersect(col_order2, colnames(tab2_data))
tab2_data <- tab2_data[, col_order2]

# Round for printing
tab2_print <- tab2_data %>% mutate(across(where(is.numeric), ~round(., 1)))
print(as.data.frame(tab2_print))
cat("\n")