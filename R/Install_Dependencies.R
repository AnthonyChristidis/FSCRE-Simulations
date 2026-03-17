# ==============================================================================
# Dependency Installation Script for FSCRE Simulations and Applications
# ==============================================================================

cat("Checking and installing required packages...\n\n")

# Set a default CRAN mirror so the script doesn't prompt the user
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 1. CRAN Packages
cran_packages <- c(
  "MASS", "Matrix", "lattice", "ggplot2", "caret", "matrixStats", "generics",
  "glmnet", "cellWise", "randomGLM", "pense", "robustbase", "robustHD", "perry",
  "doParallel", "foreach", "iterators", "remotes"
)

# Install missing CRAN packages
missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran)
} else {
  cat("All standard CRAN packages are already installed.\n")
}

# 2. Bioconductor Packages (for TCGA Application)
bioc_packages <- c(
  "BiocGenerics", "S4Vectors", "IRanges", "GenomicRanges", "Seqinfo",
  "Biobase", "SummarizedExperiment", "MultiAssayExperiment", 
  "curatedTCGAData", "TCGAutils", "MatrixGenerics"
)

# Install BiocManager if needed, then install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

missing_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages are already installed.\n")
}

# 3. GitHub Packages (CR-Lasso implementation)
if (!requireNamespace("regcell", quietly = TRUE)) {
  cat("Installing regcell from GitHub...\n")
  remotes::install_github("PengSU517/regcell")
} else {
  cat("Package 'regcell' is already installed.\n")
}

# 4. The Proposed Method (srlars)
if (!requireNamespace("srlars", quietly = TRUE)) {
  cat("Installing srlars...\n")
  install.packages("srlars")
  
  # Note: If the package is not yet approved on CRAN, use the GitHub line below instead:
  # remotes::install_github("AnthonyChristidis/srlars")
} else {
  cat("Package 'srlars' is already installed.\n")
}

cat("\n==============================================================================\n")
cat("Environment setup complete! You are ready to run the simulations.\n")
cat("==============================================================================\n")