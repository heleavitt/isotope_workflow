#!/bin/bash
#SBATCH --job-name=SCALED_mix_hypervolume
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  
#SBATCH --mem=16gb
#SBATCH --time=2-12:00:00
#SBATCH --output=scripts/logs/mix_hypervolume_%A-%a.out
#SBATCH --error=scripts/logs/mix_hypervolume_%A-%a.err
#SBATCH --array=1-1       # wrapper will set this to match #lines in combinations.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hl51981@uga.edu

cd /work/jnlab/fourchon_isotope_runs/fa22_extra_error_BA/ || {
  echo "Error: Working directory not found!" >&2
  exit 1
}

module purge
module load foss/2022a R/4.3.1-foss-2022a JAGS/4.3.1-foss-2022a

# Read the current combination line
COMBINATION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/combinations.txt)

if [ -z "$COMBINATION" ]; then
  echo "Error: No combination found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

# **Split into four variables**: species, tef, edge, buf
IFS=',' read -r SPECIES_CODE TEF_DESCRIPTOR EDGE_VAL BUF_VAL <<< "$COMBINATION"

echo "Processing species:     $SPECIES_CODE"
echo "Using TEF descriptor:   $TEF_DESCRIPTOR"
echo "Edge distance (m):      $EDGE_VAL"
echo "Buffer radius (m):      $BUF_VAL"

mkdir -p results models scripts/logs

# Pass all four arguments into your R script
Rscript scripts/fourchon_mixing_hv_scales.r \
  "$SPECIES_CODE" \
  "$TEF_DESCRIPTOR" \
  "$EDGE_VAL" \
  "$BUF_VAL" \
|| {
  echo "Error: R script execution failed!" >&2
  exit 1
}
