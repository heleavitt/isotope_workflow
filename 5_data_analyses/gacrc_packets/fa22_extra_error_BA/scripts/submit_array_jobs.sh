#!/bin/bash

# Navigate to the project directory
cd /work/jnlab/fourchon_isotope_runs/fa22_extra_error_BA/ || exit 1

# Calculate the number of lines (combinations)
ARRAY_SIZE=$(wc -l < scripts/combinations.txt)

# Check if combinations.txt is empty
if [ "$ARRAY_SIZE" -eq 0 ]; then
  echo "Error: combinations.txt is empty!" >&2
  exit 1
fi

# Submit the Slurm job with the correct array size
echo "Submitting job with array size: $ARRAY_SIZE"
sbatch --array=1-"$ARRAY_SIZE" scripts/fourchon_mixing_model_hypervolume.sh
