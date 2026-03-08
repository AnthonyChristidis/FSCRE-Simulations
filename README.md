# Fast and Scalable Cellwise-Robust Ensembles (FSCRE) - Simulations

This repository contains the R scripts and code required to reproduce the simulation studies and data analyses presented in the paper:

> **"Fast and Scalable Cellwise-Robust Ensembles for High-Dimensional Data"** 
> *Anthony-Alexander Christidis, Jeyshinee Pyneeandee, Gabriela Cohen-Freue* (Under Review)

The core methodology (the FSCRE algorithm) is implemented in the `srlars` R package, which is available on CRAN. This repository provides the scaffolding to generate the complex cellwise and casewise contamination scenarios, run the state-of-the-art competitors, and compile the performance metrics (MSPE, Precision, Recall, CPU Time).

## Repository Structure

*   **`R/`**: Contains the core simulation scripts:
    *   `generateData.R`: Generates high-dimensional data with block-collinearity and applies 5 distinct contamination scenarios (Casewise, Cellwise Marginal, Cellwise Correlation, and Mixtures).
    *   `generatePred.R`: Fits the FSCRE algorithm alongside all baseline and state-of-the-art competitor methods, returning performance metrics.
    *   `simFunc.R` & `generateOutput.R`: Wrapper functions to iterate over sparsity levels and contamination proportions.
    *   `Test_Runner.R`: A lightweight script to verify the pipeline locally.
*   **`SparseShootingS/`**: Contains the author-provided implementation of the Sparse Shooting S-estimator.
*   **`CRLasso/`**: Contains wrappers/scripts for the CR-Lasso implementation.

## Prerequisites and Installation

To run these simulations, you must install the `srlars` package along with several competitor and utility packages. 

Run the following in R to set up your environment:

```r
# 1. Install the proposed method from CRAN (or GitHub)
install.packages("srlars")
# devtools::install_github("AnthonyChristidis/srlars") # Development version

# 2. Install CRAN dependencies and baselines
install.packages(c("mvnfast", "glmnet", "cellWise", "randomGLM", "parallel", "pense", "robustHD"))

# 3. Install GitHub dependencies for competitors
# install.packages("devtools")
devtools::install_github("PengSU517/regcell") # For CR-Lasso
```

## Reproducing the Simulations

### 1. Fast Local Verification
To ensure all packages are installed correctly and the data pipelines work on your machine, you can run the localized test script. This runs a miniaturized version of the simulation ($n=30, p=50, N=2$) across all scenarios.

```bash
Rscript R/Test_Runner.R
```

This will output test results into a local `test_results/` directory.

### 2. Full Simulation Execution
The full simulation study is computationally intensive due to the high-dimensional setting ($p=500$) and the number of replications ($N=50$). To reproduce the full results, researchers should adapt the wrapper functions (`generateOutput.R`) to run in parallel on a high-performance computing environment. 

## Citation

If you use this code or the `srlars` package in your research, please cite the corresponding paper:

```bibtex
@article{christidis2026fscre,
  title={Fast and Scalable Cellwise-Robust Ensembles for High-Dimensional Data},
  author={Christidis, Anthony-Alexander and Pyneeandee, Jeyshinee and Cohen-Freue, Gabriela},
  journal={Under Review},
  year={2026}
}
```

## License

This repository is licensed under GPL (>= 2).








