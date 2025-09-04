#!/bin/bash
cd /work/jnlab/fourchon_isotope_runs/fa22_isotope_scales_extra_error/ || exit 1

SPECIES_CODES=( pensets calsap palsp )
mkdir -p scripts/logs
> scripts/combinations.txt

for species in "${SPECIES_CODES[@]}"; do
  echo "Searching TEF & MIX for: $species"

  # TEF files via glob
  tef_glob="data/tef_files/tef_${species}_step_*.csv"
  shopt -s nullglob
  tef_list=( $tef_glob )
  shopt -u nullglob

  if [ ${#tef_list[@]} -eq 0 ]; then
    echo "Warning: no TEF for $species" >> scripts/logs/missing_files.log
    continue
  fi

  # MIX files via glob
  mix_glob="data/mix_files/edge_*_buf_*_${species}_mix.csv"
  shopt -s nullglob
  mix_list=( $mix_glob )
  shopt -u nullglob

  if [ ${#mix_list[@]} -eq 0 ]; then
    echo "Warning: no MIX for $species" >> scripts/logs/missing_files.log
    continue
  fi

  # Nested loops
  for tef_file in "${tef_list[@]}"; do
    tef_desc=$(basename "$tef_file" \
                 | sed -e "s/tef_${species}_step_//" -e 's/\.csv//')
    echo "  TEF step: $tef_desc"

    for mix_file in "${mix_list[@]}"; do
      mix_base=$(basename "$mix_file")
      edge_val=$(echo "$mix_base" \
                   | sed -n 's/edge_\([0-9]*\)_buf.*/\1/p')
      buf_val=$(echo "$mix_base" \
                  | sed -n 's/.*_buf_\([0-9]*\)_'"${species}"'_mix\.csv/\1/p')

      if [[ -n $edge_val && -n $buf_val ]]; then
        echo "${species},${tef_desc},${edge_val},${buf_val}" \
          >> scripts/combinations.txt
      else
        echo "Warn: parse failed for $mix_base" \
          >> scripts/logs/missing_files.log
      fi
    done
  done
done

echo "✅ Done — see scripts/combinations.txt"
