#!/bin/bash
#SBATCH --job-name=tc-grid-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --clusters=all
#SBATCH --partition=long
#SBATCH --output=jobs/log/phase0_%j.out
#SBATCH --error=jobs/log/phase0_%j.err

# Phase 0: Preprocess
# Generate networks, rasterised grids, and downscaling factors for all nations
# Targets:
# - results/power/by_country/<ISO_A3>/exposure/edges_split.geoparquet
# - results/power/by_country/<ISO_A3>/storms/downscale_factors.npy

set -euo pipefail

echo "Starting Phase 0: Preprocess"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"

# Ensure log directory exists
mkdir -p jobs/log

# Generate target list dynamically
TARGETS_FILE="jobs/log/phase0_targets_${SLURM_JOB_ID}.txt"
echo "Generating target list..."

# Read countries from config
COUNTRIES=($(jq -r '.[]' config/tc_grid/permitted_countries.json))

# Generate targets for each country
> "$TARGETS_FILE"  # Clear/create file
for COUNTRY in "${COUNTRIES[@]}"; do
    echo "results/power/by_country/${COUNTRY}/exposure/edges_split.geoparquet" >> "$TARGETS_FILE"
    echo "results/power/by_country/${COUNTRY}/storms/downscale_factors.npy" >> "$TARGETS_FILE"
done

echo "Generated $(wc -l < "$TARGETS_FILE") targets"

# Add sliced track targets for all countries, storm sets, and samples
echo "Adding sliced track targets..."
STORM_SETS=($(jq -r '.[]' config/tc_grid/storm_sets.json))
for COUNTRY in "${COUNTRIES[@]}"; do
    for STORM_SET in "${STORM_SETS[@]}"; do
        for SAMPLE in 0 1 2 3 4; do
            echo "results/power/by_country/${COUNTRY}/storms/${STORM_SET}/${SAMPLE}/tracks.geoparquet" >> "$TARGETS_FILE"
        done
    done
done

echo "Total targets including sliced tracks: $(wc -l < "$TARGETS_FILE")"

# Unlock snakemake directory
pixi run snakemake --cores 1 --unlock

# Run snakemake with whole node resources
# Build all targets without rule restrictions
# Phase 1 will be responsible for only running estimate_wind_fields and electricity_grid_damages
pixi run snakemake \
    --cores 40 \
    --resources mem_mb=140000 \
    --rerun-incomplete \
    $(cat "$TARGETS_FILE")

EXIT_CODE=$?

# Cleanup targets file
rm -f "$TARGETS_FILE"

echo "Phase 0 completed at: $(date)"
exit $EXIT_CODE
