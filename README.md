# Fast and Scalable Cellwise-Robust Ensembles (FSCRE) - Simulations

This repository contains the R scripts and code required to reproduce the simulation studies and data analyses presented in the paper:

> **"Fast and Scalable Cellwise-Robust Ensembles for High-Dimensional Data"** 
> *Anthony-Alexander Christidis, Jeyshinee Pyneeandee, Gabriela Cohen-Freue* (Under Review)

The core methodology (the FSCRE algorithm) is implemented in the `srlars` R package, which is available on CRAN. This repository provides the scaffolding to generate the complex cellwise and casewise contamination scenarios, run the state-of-the-art competitors, and compile the performance metrics (MSPE, Precision, Recall, CPU Time) for both synthetic simulations and real-world bioinformatics applications.

## Repository Structure

*   **`R/`**: Contains the core simulation and application scripts:
    *   `generateData.R`: Generates high-dimensional data with block-collinearity and applies 5 distinct contamination scenarios (Casewise, Cellwise Marginal, Cellwise Correlation, and Mixtures).
    *   `generatePred.R`: Fits the FSCRE algorithm alongside all baseline and state-of-the-art competitor methods, returning performance metrics.
    *   `simFunc.R` & `generateOutput.R`: Wrapper functions to iterate over sparsity levels and contamination proportions.
    *   `Generate_Results.R`: The main execution script for the full simulation study. It iterates over all contamination scenarios, SNRs, and parameter grids to generate the complete set of results.
    *   `Test_Runner.R`: A lightweight script to verify the simulation pipeline locally.
    *   `Generate_CPU_Results.R`: Script to conduct the computational scalability study, measuring execution time across varying dimensions ($p$) and sample sizes ($n$).
    *   `Application_TCGA.R`: Script to reproduce the Proteogenomics (TCGA BRCA) real data application (predicting protein abundance from mRNA).
*   **`SparseShootingS/`**: Contains the author-provided implementation of the Sparse Shooting S-estimator.

## Prerequisites and Installation

To ensure a smooth setup, we provide a single script that automatically checks for and installs all required dependencies from CRAN, Bioconductor, and GitHub (including the `srlars` package and competitor implementations).

Run the following command in your R console from the root directory of this repository:

```r
source("Install_Dependencies.R")
```

## Reproducing the Analyses

### 1. Fast Local Verification (Simulations)
To ensure all packages are installed correctly and the data pipelines work on your machine, you can run the localized test script. This runs a miniaturized version of the simulation ($n=30, p=50, N=2$) across all scenarios.

```bash
Rscript R/Test_Runner.R
```

This will output test results into a local `test_results/` directory.

### 2. Full Simulation Execution
The full simulation study is computationally intensive due to the high-dimensional setting ($p=500$) and the number of replications ($N=50$). To reproduce the full results, researchers should adapt the wrapper functions (`generateOutput.R`) to run in parallel on a high-performance computing environment. 

### 3. Computational Scalability Study
To reproduce the CPU timing results demonstrating the scalability of FSCRE against increasing dimensions and sample sizes, run the timing script. This isolates the most challenging contamination scenario to stress-test the algorithms.

```bash
Rscript R/Generate_CPU_Results.R
```

### 4. Bioinformatics Data Application
To reproduce the empirical results on real-world genomic data, run the application script. This script automatically downloads the necessary multi-omics data, performs intersections, introduces targeted artificial contamination for robustness testing, and evaluates the models.

```bash
# Run Proteogenomics application
Rscript R/Application_TCGA.R
```

Results will be saved as `.rds` files in the `results/` directory.

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