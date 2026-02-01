#! /bin/bash -l
#SBATCH --job-name="New Sim Study"
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -t 330:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lapo.santi@ucdconnect.ie

# Always good hygiene on clusters
module purge

# 1) Load GCC toolchain (g++ 13.2.0 + runtime)
module load gcc/13.2.0-gcc-11.5.0-yrgytpa

# 2) Load R AFTER the compiler so it sees the same toolchain
module load R/4.4.2

# 3) (Optional but recommended) load git
module load git || true

# 4) Make sure we use the matching libstdc++ at runtime
LIBSTDCPP_13=$(g++ -print-file-name=libstdc++.so.6)
export LD_LIBRARY_PATH="$(dirname "$LIBSTDCPP_13"):$LD_LIBRARY_PATH"

cd "$SLURM_SUBMIT_DIR"

echo "[$(date)] Pulling latest main..."
git pull origin main

echo "[$(date)] Starting R script..."
Rscript new_sim_study.R
echo "[$(date)] Done."
