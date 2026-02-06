#! /bin/bash -l
#SBATCH --job-name="PE Comparison"
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lapo.santi@ucdconnect.ie
#SBATCH --output=logs/pe_comparison_%j.out
#SBATCH --error=logs/pe_comparison_%j.err

set -euo pipefail

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
export LD_LIBRARY_PATH="$(dirname "$LIBSTDCPP_13"):${LD_LIBRARY_PATH:-}"

# Avoid thread oversubscription
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR:-$PWD}"

mkdir -p logs results

echo "[$(date)] Pulling latest main..."
git pull origin main || true

# Defaults match Comparison with Pearce_Ereshova.R (and new_sim_study.R seed base)
export SEED="${SEED:-123}"
export GAMMA_TRUE="${GAMMA_TRUE:-0.71}"
export P_ADJ="${P_ADJ:-0.85}"
export SIM_DESIGN="${SIM_DESIGN:-}"
export SIM_K_TRUE="${SIM_K_TRUE:-3}"

export T_ITER="${T_ITER:-3000}"
export T_BURN="${T_BURN:-600}"
export N_CHAINS="${N_CHAINS:-4}"
export USE_PARALLEL="${USE_PARALLEL:-TRUE}"

log_tag="design${SIM_DESIGN:-default}_K${SIM_K_TRUE}_seed${SEED}_job${SLURM_JOB_ID:-local}"

echo "[$(date)] Starting Comparison with Pearce_Ereshova.R (${log_tag})..."
Rscript "Comparison with Pearce_Ereshova.R" 2>&1 | tee "logs/pe_comparison_${log_tag}.log"
echo "[$(date)] Done."
