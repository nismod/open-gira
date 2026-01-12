#!/usr/bin/env python3
"""
Generate phase1 job scripts for small, medium, and large countries.
Dynamically determines country categories based on target counts.
"""

import argparse
import json
import csv
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Category thresholds (number of gridfinder targets for each country)
SMALL_THRESHOLD = 4000
LARGE_THRESHOLD = 10000

# Resource specifications
RESOURCES = {
    "small": {"cpus": 2, "mem": "6G", "time": "2:00:00", "max_concurrent": 100},
    "medium": {"cpus": 8, "mem": "24G", "time": "4:00:00", "max_concurrent": 50},
    "large": {"cpus": 40, "mem": "120G", "time": "6:00:00", "max_concurrent": 20},
}

SCRIPT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=tc-grid-phase1-{category}
#SBATCH --array=0-{max_array_id}%{max_concurrent}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --partition=long
#SBATCH --output=jobs/log/phase1_{category}_%A_%a.out
#SBATCH --error=jobs/log/phase1_{category}_%A_%a.err

# Phase 1 {category_title}: Core processing for {category} countries (target_count {threshold_desc})
# Each array job handles one (country, storm_set, sample) combination
# Resources: {cpus} CPU, {mem} memory
# Targets:
# - results/power/by_country/<ISO_A3>/disruption/<ATMOSPHERE>/<SAMPLE>

set -euo pipefail

echo "Starting Phase 1 {category_title}: Core processing"
echo "Job ID: $SLURM_JOB_ID"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"

# Ensure log directory exists
mkdir -p jobs/log

# Base directory (where the original workspace is)
BASE_DIR="/data/ouce-opsis/cenv0899/2025-tc-grid/open-gira"

# {category_title} countries (target_count {threshold_desc}) - permitted countries only
COUNTRIES=({countries})

STORM_SETS=($(jq -r '.[]' "${{BASE_DIR}}/config/tc_grid/storm_sets.json"))
SAMPLES=({samples})

NUM_COUNTRIES=${{#COUNTRIES[@]}}
NUM_STORM_SETS=${{#STORM_SETS[@]}}
NUM_SAMPLES=${{#SAMPLES[@]}}

# Calculate total combinations
TOTAL_COMBINATIONS=$((NUM_COUNTRIES * NUM_STORM_SETS * NUM_SAMPLES))

echo "Total combinations ({category}): $TOTAL_COMBINATIONS"
echo "Countries: $NUM_COUNTRIES, Storm sets: $NUM_STORM_SETS, Samples: $NUM_SAMPLES"

# Map array task ID to (country, storm_set, sample)
TASK_ID=$SLURM_ARRAY_TASK_ID

COUNTRY_IDX=$((TASK_ID / (NUM_STORM_SETS * NUM_SAMPLES)))
REMAINDER=$((TASK_ID % (NUM_STORM_SETS * NUM_SAMPLES)))
STORM_SET_IDX=$((REMAINDER / NUM_SAMPLES))
SAMPLE_IDX=$((REMAINDER % NUM_SAMPLES))

COUNTRY=${{COUNTRIES[$COUNTRY_IDX]}}
STORM_SET=${{STORM_SETS[$STORM_SET_IDX]}}
SAMPLE=${{SAMPLES[$SAMPLE_IDX]}}

echo "Processing: Country=$COUNTRY, Storm_set=$STORM_SET, Sample=$SAMPLE"

# Create unique working directory for this task
WORK_DIR="${{BASE_DIR}}_work_{category}_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"
echo "Working directory: $WORK_DIR"

# Copy workspace (excluding heavy/shared directories)
# .pixi is excluded to avoid copying GBs of Python environments (use shared environment)
echo "Copying workspace to working directory..."
mkdir -p "$WORK_DIR"
rsync -a \\
    --exclude='.snakemake' \\
    --exclude='.git' \\
    --exclude='.pixi' \\
    --exclude='results' \\
    --exclude='open-gira_work_*' \\
    "${{BASE_DIR}}/" "$WORK_DIR/"

# Link to shared directories
ln -sf "${{BASE_DIR}}/results" "$WORK_DIR/results"
ln -sf "${{BASE_DIR}}/.pixi" "$WORK_DIR/.pixi"

# Change to working directory
cd "$WORK_DIR"

# Define the target for this specific combination
TARGET="results/power/by_country/${{COUNTRY}}/disruption/${{STORM_SET}}/${{SAMPLE}}"

echo "Target: $TARGET"

# Unlock snakemake directory
pixi run snakemake --cores 1 --unlock

# Run snakemake with strict rule limiting to Phase 1 rules only
# All preprocessing (networks, downscaling, tracks) must be done in Phase 0
# This ensures embarrassingly parallel execution with no shared file creation
pixi run snakemake \\
    --allowed-rules estimate_wind_fields electricity_grid_damages \\
    --cores {cpus} \\
    --rerun-incomplete \\
    "$TARGET"

EXIT_CODE=$?

# Cleanup working directory
echo "Cleaning up working directory..."
cd "${{BASE_DIR}}"
# Remove symlinks explicitly (extra safety to ensure we don't follow them)
rm -f "$WORK_DIR/results"
rm -f "$WORK_DIR/.pixi"
# Remove the working directory
rm -rf "$WORK_DIR"

echo "Phase 1 {category_title} task completed at: $(date)"
exit $EXIT_CODE
"""


def load_target_counts(csv_path):
    """Load target counts from CSV file."""
    target_counts = {}
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            iso_a3 = row["iso_a3"]
            count = int(row["target_count"])
            target_counts[iso_a3] = count
    return target_counts


def load_permitted_countries(json_path):
    """Load permitted countries from JSON file."""
    with open(json_path, "r") as f:
        return set(json.load(f))


def load_storm_sets(json_path):
    """Load storm sets from JSON file."""
    with open(json_path, "r") as f:
        return json.load(f)


def categorize_countries(target_counts, permitted_countries):
    """Categorize countries into small, medium, and large."""
    categories = {"small": [], "medium": [], "large": []}

    for country in permitted_countries:
        if country not in target_counts:
            logger.warning(
                f"{country} in permitted list but not in target counts, skipping"
            )
            continue

        count = target_counts[country]
        if count < SMALL_THRESHOLD:
            categories["small"].append(country)
        elif count < LARGE_THRESHOLD:
            categories["medium"].append(country)
        else:
            categories["large"].append(country)

    # Sort for consistency
    for cat in categories:
        categories[cat].sort()

    return categories


def generate_script(category, countries, samples, num_storm_sets, output_path):
    """Generate a phase 1 script for the given category."""
    resources = RESOURCES[category]

    # Calculate array size
    num_countries = len(countries)
    num_samples = len(samples)
    total_combinations = num_countries * num_storm_sets * num_samples
    max_array_id = total_combinations - 1

    # Determine threshold description
    if category == "small":
        threshold_desc = f"< {SMALL_THRESHOLD}"
    elif category == "medium":
        threshold_desc = f">= {SMALL_THRESHOLD} and < {LARGE_THRESHOLD}"
    else:  # large
        threshold_desc = f">= {LARGE_THRESHOLD}"

    # Format countries and samples for bash array
    countries_str = " ".join(countries)
    samples_str = " ".join(map(str, samples))

    # Generate script from template
    script_content = SCRIPT_TEMPLATE.format(
        category=category,
        category_title=category.capitalize(),
        cpus=resources["cpus"],
        mem=resources["mem"],
        time=resources["time"],
        max_concurrent=resources["max_concurrent"],
        max_array_id=max_array_id,
        threshold_desc=threshold_desc,
        countries=countries_str,
        samples=samples_str,
    )

    # Write script
    output_path.write_text(script_content)
    output_path.chmod(0o755)

    logger.info(f"Generated {output_path.name}:")
    logger.info(f"  Countries: {num_countries}")
    logger.info(f"  Samples: {num_samples}")
    logger.info(f"  Storm sets: {num_storm_sets}")
    logger.info(f"  Total array tasks: {total_combinations}")
    logger.info(f"  Resources: {resources['cpus']} CPU, {resources['mem']} mem")
    logger.info("")


def main():
    parser = argparse.ArgumentParser(
        description="Generate phase1 job scripts for TC-grid workflow"
    )
    parser.add_argument(
        "--samples",
        type=str,
        default="0,1,2,3,4",
        help="Comma-separated list of sample numbers (default: 0,1,2,3,4)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("./jobs/tc-grid/"),
        help="Output directory for generated scripts (default: ./jobs/tc-grid/)",
    )

    args = parser.parse_args()

    # Parse samples
    samples = [int(s.strip()) for s in args.samples.split(",")]
    logger.info(f"Samples to process: {samples}")
    logger.info("")

    # Load data
    base_dir = Path(__file__).parent.parent.parent
    target_counts = load_target_counts(
        base_dir / "results/power/target_count_by_country.csv"
    )
    permitted_countries = load_permitted_countries(
        base_dir / "config/tc_grid/permitted_countries.json"
    )
    storm_sets = load_storm_sets(base_dir / "config/tc_grid/storm_sets.json")

    num_storm_sets = len(storm_sets)
    logger.info(f"Loaded {len(target_counts)} countries with target counts")
    logger.info(f"Loaded {len(permitted_countries)} permitted countries")
    logger.info(f"Loaded {num_storm_sets} storm sets")
    logger.info("")

    # Categorize countries
    categories = categorize_countries(target_counts, permitted_countries)

    logger.info("Country categorization:")
    logger.info(f"  Small (< {SMALL_THRESHOLD}): {len(categories['small'])} countries")
    logger.info(
        f"  Medium ({SMALL_THRESHOLD}-{LARGE_THRESHOLD}): {len(categories['medium'])} countries"
    )
    logger.info(f"  Large (>= {LARGE_THRESHOLD}): {len(categories['large'])} countries")
    logger.info("")

    # Generate scripts
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for category in ["small", "medium", "large"]:
        if not categories[category]:
            logger.warning(
                f"No countries in {category} category, skipping script generation"
            )
            continue

        output_path = output_dir / f"1-core-{category}.sh"
        generate_script(
            category, categories[category], samples, num_storm_sets, output_path
        )

    logger.info("Script generation complete!")


if __name__ == "__main__":
    main()
