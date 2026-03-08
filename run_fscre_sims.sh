#!/bin/bash
#SBATCH --job-name=fscre_simulation                     # Job name
#SBATCH --output=sim_study_out/sim_study_%A_%a.out      # Name of stdout output file (%A = array job ID, %a = task ID)
#SBATCH --error=sim_study_err/sim_study_%A_%a.err       # Name of stderr error file
#SBATCH --partition=short                               # Short partition
#SBATCH --time=6:00:00                                  # 12 hours should be plenty for one scenario/SNR combo
#SBATCH --ntasks=1                                      # Number of tasks (always 1 for R scripts)
#SBATCH --cpus-per-task=6                               # Matches makeCluster(5) in generatePred.R
#SBATCH --mem=12G                                       # Ample RAM for job
#SBATCH --array=1-15                                    # Array job for 15 combinations (5 scenarios * 3 SNRs)
#SBATCH --mail-type=END,FAIL                            # Send email when job ends or fails
#SBATCH --mail-user=anthony-alexander_christidis@hms.harvard.edu
  
# Print job information
echo "Starting job at: $(date)"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Working directory: $SLURM_SUBMIT_DIR"
echo "Available cores: $SLURM_CPUS_PER_TASK"
echo "Requested memory: ${SLURM_MEM_PER_NODE}MB"

# Ensure we're in the submission directory
cd $SLURM_SUBMIT_DIR

# Create necessary directories if they don't exist
mkdir -p results
mkdir -p sim_study_out
mkdir -p sim_study_err

# Set thread limits for numerical libraries
# CRITICAL: Prevents openBLAS/MKL from trying to use all node cores,
# which conflicts with our explicit makeCluster(5) parallelism.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export LAPACK_NUM_THREADS=1

# Note: Adjust R_LIBS_USER if your packages are installed elsewhere for R 4.4.2
# export R_LIBS_USER="$HOME/R/4.4.2/library"

# Load the gcc and R modules (Verify these exist on O2 currently)
module load gcc/14.2.0 R/4.4.2

# Print system resources
echo "System info:"
free -h | head -2
echo ""

# Run the simulation
# Executing the array runner script located in the R/ folder
echo "Starting R script at: $(date)"
Rscript R/Generate_Results_Array.R

echo "Job completed at: $(date)"

# Print final file info
echo "Results should be in:"
ls -la results/ | tail -5