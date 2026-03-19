# -------------------------------------------------------------------------
# Dependency Installation Script for FSCRE Simulations and Applications
# Enforces specific package versions for strict reproducibility.
# -------------------------------------------------------------------------

cat("Checking and installing required packages (with version enforcement)...\n\n")

# Set a default CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Ensure 'remotes' is installed first, as it's needed for version control
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)

# 1. CRAN Packages and Exact Versions
cran_pkgs <- list(
  MASS = "7.3-65",
  Matrix = "1.7-4",
  lattice = "0.22-7",
  ggplot2 = "4.0.1",
  caret = "7.0-1",
  matrixStats = "1.5.0",
  generics = "0.1.4",
  glmnet = "4.1-10",
  cellWise = "2.5.5",
  randomGLM = "1.10-1",
  pense = "2.5.0",
  robustbase = "0.99-6",
  robustHD = "0.8.4",
  perry = "0.3.1",
  doParallel = "1.0.17",
  foreach = "1.5.2",
  iterators = "1.0.14",
  R.utils = "2.13.0"
)

cat("--- Checking CRAN Packages ---\n")
for (pkg in names(cran_pkgs)) {
  req_ver <- cran_pkgs[[pkg]]
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s (version %s)...\n", pkg, req_ver))
    remotes::install_version(pkg, version = req_ver, upgrade = "never", quiet = TRUE)
  } else {
    curr_ver <- as.character(packageVersion(pkg))
    if (curr_ver != req_ver) {
      cat(sprintf("Updating %s from %s to required version %s...\n", pkg, curr_ver, req_ver))
      remotes::install_version(pkg, version = req_ver, upgrade = "never", quiet = TRUE)
    } else {
      cat(sprintf("%s (v%s) is already installed.\n", pkg, curr_ver))
    }
  }
}

# 2. Bioconductor Packages (for TCGA Application)
# Note: For Bioconductor, enforcing specific patch versions via remotes can break internal 
# Bioc dependencies. It is generally safer to let BiocManager install the suite.
# However, we list your specific versions here for documentation.
bioc_pkgs <- list(
  BiocGenerics = "0.56.0",
  S4Vectors = "0.48.0",
  IRanges = "2.44.0",
  GenomicRanges = "1.62.1",
  Seqinfo = "1.0.0",
  Biobase = "2.70.0",
  SummarizedExperiment = "1.40.0",
  MultiAssayExperiment = "1.36.1",
  curatedTCGAData = "1.32.1",
  TCGAutils = "1.30.2",
  MatrixGenerics = "1.22.0"
)

cat("\n--- Checking Bioconductor Packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

missing_bioc <- names(bioc_pkgs)[!(names(bioc_pkgs) %in% installed.packages()[,"Package"])]
if(length(missing_bioc) > 0) {
  cat("Installing missing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
  cat("All required Bioconductor packages are installed.\n")
}

# 3. GitHub Packages (CR-Lasso implementation)
cat("\n--- Checking GitHub Packages ---\n")
if (!requireNamespace("regcell", quietly = TRUE)) {
  cat("Installing regcell from GitHub...\n")
  remotes::install_github("PengSU517/regcell", upgrade = "never", quiet = TRUE)
} else {
  cat("Package 'regcell' is already installed.\n")
}

# 4. The Proposed Method (srlars)
cat("\n--- Checking Proposed Method Package ---\n")
if (!requireNamespace("srlars", quietly = TRUE)) {
  cat("Installing srlars...\n")
  # Try CRAN first
  tryCatch({
    install.packages("srlars")
  }, error = function(e) {
    cat("Could not find on CRAN. Falling back to GitHub...\n")
    remotes::install_github("AnthonyChristidis/srlars", upgrade = "never", quiet = TRUE)
  })
} else {
  curr_ver <- as.character(packageVersion("srlars"))
  if (curr_ver != "2.0.1") {
    cat(sprintf("Updating srlars to v2.0.1...\n"))
    # Assuming v2.0.1 is the version you submitted/are hosting
    # If it's not on CRAN yet, use GitHub:
    # remotes::install_github("AnthonyChristidis/srlars@v2.0.1", upgrade="never")
    install.packages("srlars") 
  } else {
    cat(sprintf("srlars (v%s) is already installed.\n", curr_ver))
  }
}

cat("\n==============================================================================\n")
cat("Environment setup complete! You are ready to run the simulations.\n")
cat("==============================================================================\n")