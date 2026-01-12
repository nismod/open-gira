#!/bin/bash
#SBATCH --job-name=tc-grid-phase2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --partition=long
#SBATCH --output=jobs/log/phase2_%j.out
#SBATCH --error=jobs/log/phase2_%j.err

# Phase 2: Postprocess
# Aggregate and analyse results
# Targets:
# - results/power/by_storm_set/<ATMOSPHERE>/disruption/EAPA_admin-level-2.gpq
# - results/power/by_storm_set/<ATMOSPHERE>/disruption/pop_affected_RP

set -euo pipefail

echo "Starting Phase 2: Postprocess"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"

# Ensure log directory exists
mkdir -p jobs/log

# Generate target list dynamically
TARGETS_FILE="jobs/log/phase2_targets_${SLURM_JOB_ID}.txt"
echo "Generating target list..."

# Read storm sets from config
STORM_SETS=($(jq -r '.[]' config/tc_grid/storm_sets.json))

# Generate targets for each storm set
> "$TARGETS_FILE"  # Clear/create file
for STORM_SET in "${STORM_SETS[@]}"; do
    echo "results/power/by_storm_set/${STORM_SET}/disruption/EAPA_admin-level-2.gpq" >> "$TARGETS_FILE"
    echo "results/power/by_storm_set/${STORM_SET}/disruption/pop_affected_RP" >> "$TARGETS_FILE"
done

echo "Generated $(wc -l < "$TARGETS_FILE") targets"

# Unlock snakemake directory
pixi run snakemake --cores 1 --unlock

# Run snakemake with whole node resources
# Strictly limit to post-processing rules (everything except Phase 1 rules)
# We use --omit-from to exclude the Phase 1 rules
pixi run snakemake \
    --cores 40 \
    --resources mem_mb=140000 \
    --omit-from estimate_wind_fields \
    --omit-from electricity_grid_damages \
    --rerun-incomplete \
    $(cat "$TARGETS_FILE")

EXIT_CODE=$?

# Cleanup targets file
rm -f "$TARGETS_FILE"

echo "Phase 2 completed at: $(date)"
exit $EXIT_CODE
