#!/usr/bin/env python3
"""
Generate phase1 job scripts for small, medium, and large countries.
Dynamically determines country categories based on target counts.
"""

import argparse
import json
import csv
import logging
import subprocess
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Category thresholds (number of gridfinder targets for each country)
SMALL_THRESHOLD = 4000
LARGE_THRESHOLD = 30000

# Resource specifications
RESOURCES = {
    "small": {"cpus": 2, "mem": "6G", "time": "4:00:00", "max_concurrent": 100},
    "medium": {"cpus": 8, "mem": "60G", "time": "16:00:00", "max_concurrent": 50},
    "large": {"cpus": 96, "mem": "230G", "time": "48:00:00", "max_concurrent": 20},
}

SCRIPT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=tc-grid-1-{category}
#SBATCH --array=0-{max_array_id}%{max_concurrent}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --clusters=all
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

{array_setup}

# Map array task ID to (country, storm_set, sample)
TASK_ID=$SLURM_ARRAY_TASK_ID

{index_calculation}

COUNTRY={country_lookup}
STORM_SET={storm_set_lookup}
SAMPLE={sample_lookup}

echo "Processing: Country=$COUNTRY, Storm_set=$STORM_SET, Sample=$SAMPLE"

# Change to base directory
cd "$BASE_DIR"

# Capture path to snakemake
SNAKEMAKE=$(pixi run which snakemake)

# Define the target for this specific combination
TARGET="results/power/by_country/${{COUNTRY}}/disruption/${{STORM_SET}}/${{SAMPLE}}"

echo "Target: $TARGET"

# Run snakemake with strict rule limiting to Phase 1 rules only
# All preprocessing (networks, downscaling, tracks) must be done in Phase 0
# This ensures embarrassingly parallel execution with no shared file creation
# Use --nolock to allow multiple snakemake processes in parallel
"$SNAKEMAKE" \\
    --allowed-rules estimate_wind_fields electricity_grid_damages \\
    --cores {cpus} \\
    --rerun-incomplete \\
    --nolock \\
    "$TARGET"

EXIT_CODE=$?

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


def check_missing_targets(countries, storm_sets, samples, base_dir):
    """
    Check which targets are missing or incomplete using snakemake.

    Two-phase approach:
    1. Quick filesystem check - if directory doesn't exist, it's definitely missing
    2. Snakemake check - for existing directories, verify they're complete

    Returns list of (index, country, storm_set, sample) tuples for missing targets.
    """
    logger.info("Checking for missing targets...")

    # Phase 1: Quick filesystem check
    missing_targets = []
    needs_snakemake_check = []
    target_metadata = {}  # maps target path to (index, country, storm_set, sample)

    total_targets = 0
    for country_idx, country in enumerate(countries):
        for storm_set_idx, storm_set in enumerate(storm_sets):
            for sample_idx, sample in enumerate(samples):
                index = (
                    country_idx * (len(storm_sets) * len(samples))
                    + storm_set_idx * len(samples)
                    + sample_idx
                )
                target = f"results/power/by_country/{country}/disruption/{storm_set}/{sample}"
                target_path = base_dir / target
                total_targets += 1

                # Quick check: does the directory exist?
                if not target_path.exists():
                    # Definitely missing, no need to check with snakemake
                    missing_targets.append((index, country, storm_set, sample))
                else:
                    # Directory exists, but might be incomplete - check with snakemake
                    needs_snakemake_check.append(target)
                    target_metadata[target] = (index, country, storm_set, sample)

    logger.info(f"Total target combinations: {total_targets}")
    logger.info(f"  Definitely missing (no directory): {len(missing_targets)}")
    logger.info(
        f"  Need snakemake check (directory exists): {len(needs_snakemake_check)}"
    )

    # Phase 2: Snakemake check for existing directories
    if not needs_snakemake_check:
        logger.info("All targets are missing, no snakemake checks needed")
        return missing_targets

    # Optimization: if ALL targets exist on filesystem, assume the run is complete
    if len(missing_targets) == 0:
        logger.info(
            "All target directories exist - assuming run is complete, skipping snakemake verification"
        )
        return []

    logger.info("Running snakemake to check existing targets for completeness...")

    # Check targets in batches to avoid overwhelming snakemake with huge DAG
    batch_size = 50

    for batch_start in range(0, len(needs_snakemake_check), batch_size):
        batch_end = min(batch_start + batch_size, len(needs_snakemake_check))
        batch = needs_snakemake_check[batch_start:batch_end]

        logger.info(
            f"  Batch {batch_start // batch_size + 1}/{(len(needs_snakemake_check) + batch_size - 1) // batch_size} ({len(batch)} targets)..."
        )

        # Run snakemake to check target status
        cmd = [
            "pixi",
            "run",
            "snakemake",
            "--cores",
            "1",
            "--nolock",
            "--dry-run",
            "--rerun-incomplete",
            "--summary",
        ] + batch

        try:
            result = subprocess.run(
                cmd,
                cwd=base_dir,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout per batch
            )
        except subprocess.TimeoutExpired:
            logger.error(
                f"Snakemake check timed out for batch {batch_start // batch_size + 1}"
            )
            raise

        if result.returncode != 0:
            logger.error("Snakemake check failed:")
            logger.error(result.stderr)
            raise RuntimeError("Failed to check targets with snakemake")

        # Parse the summary output
        # Format: output_file\tdate\trule\tlog-file(s)\tstatus\tplan
        # We want targets where status == "missing" or plan == "update pending"
        lines = result.stdout.strip().split("\n")
        for line in lines:
            # Skip header and non-data lines
            if line.startswith("output_file") or not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                continue

            target_path = parts[0]
            status = parts[4] if len(parts) > 4 else ""
            plan = parts[5] if len(parts) > 5 else ""

            # Check if this is one of our disruption targets
            if target_path in target_metadata:
                # Target needs building if it's missing or has update pending
                if status == "missing" or "update pending" in plan:
                    missing_targets.append(target_metadata[target_path])

    logger.info(
        f"Found {len(missing_targets)} missing/incomplete targets out of {total_targets}"
    )
    if total_targets > 0:
        logger.info(
            f"Reduction: {len(missing_targets)} / {total_targets} = {len(missing_targets) / total_targets:.1%}"
        )

    return missing_targets


def generate_script(
    category, countries, samples, storm_sets, output_path, check_missing=True
):
    """Generate a phase 1 script for the given category."""
    base_dir = Path(__file__).parent.parent.parent
    resources = RESOURCES[category]

    # Check which targets are missing
    if check_missing:
        missing_targets = check_missing_targets(
            countries, storm_sets, samples, base_dir
        )

        if not missing_targets:
            logger.info(
                f"No missing targets for {category} category, skipping script generation"
            )
            return

        # Generate mapping arrays for the missing targets only
        countries_map = " ".join(f"[{i}]={t[1]}" for i, t in enumerate(missing_targets))
        storm_sets_map = " ".join(
            f"[{i}]={t[2]}" for i, t in enumerate(missing_targets)
        )
        samples_map = " ".join(f"[{i}]={t[3]}" for i, t in enumerate(missing_targets))

        total_combinations = len(missing_targets)
        max_array_id = total_combinations - 1

        # Build template components for mapping mode
        array_setup = f"""# Sparse array mode: only processing missing targets
# Total missing targets: {total_combinations}

# Declare associative arrays mapping task ID to actual values
declare -A COUNTRY_MAP=({countries_map})
declare -A STORM_SET_MAP=({storm_sets_map})
declare -A SAMPLE_MAP=({samples_map})

echo "Total missing targets ({category}): {total_combinations}\""""

        index_calculation = ""  # No calculation needed
        country_lookup = "${COUNTRY_MAP[$TASK_ID]}"
        storm_set_lookup = "${STORM_SET_MAP[$TASK_ID]}"
        sample_lookup = "${SAMPLE_MAP[$TASK_ID]}"

    else:
        # Original behavior: generate all combinations
        num_countries = len(countries)
        num_storm_sets = len(storm_sets)
        num_samples = len(samples)
        total_combinations = num_countries * num_storm_sets * num_samples
        max_array_id = total_combinations - 1

        # For non-mapping mode, use original arrays
        countries_str = " ".join(countries)
        samples_str = " ".join(map(str, samples))

        # Build template components for full mode
        array_setup = f"""# {category.capitalize()} countries (target_count {{threshold_desc}}) - permitted countries only
COUNTRIES=({countries_str})

STORM_SETS=($(jq -r '.[]' "${{{{BASE_DIR}}}}/config/tc_grid/storm_sets.json"))
SAMPLES=({samples_str})

NUM_COUNTRIES=${{{{#COUNTRIES[@]}}}}
NUM_STORM_SETS=${{{{#STORM_SETS[@]}}}}
NUM_SAMPLES=${{{{#SAMPLES[@]}}}}

# Calculate total combinations
TOTAL_COMBINATIONS=$((NUM_COUNTRIES * NUM_STORM_SETS * NUM_SAMPLES))

echo "Total combinations ({category}): $TOTAL_COMBINATIONS"
echo "Countries: $NUM_COUNTRIES, Storm sets: $NUM_STORM_SETS, Samples: $NUM_SAMPLES\""""

        index_calculation = f"""COUNTRY_IDX=$((TASK_ID / (NUM_STORM_SETS * NUM_SAMPLES)))
REMAINDER=$((TASK_ID % (NUM_STORM_SETS * NUM_SAMPLES)))
STORM_SET_IDX=$((REMAINDER / NUM_SAMPLES))
SAMPLE_IDX=$((REMAINDER % NUM_SAMPLES))"""

        country_lookup = "${COUNTRIES[$COUNTRY_IDX]}"
        storm_set_lookup = "${STORM_SETS[$STORM_SET_IDX]}"
        sample_lookup = "${SAMPLES[$SAMPLE_IDX]}"

    # Determine threshold description
    if category == "small":
        threshold_desc = f"< {SMALL_THRESHOLD}"
    elif category == "medium":
        threshold_desc = f">= {SMALL_THRESHOLD} and < {LARGE_THRESHOLD}"
    else:  # large
        threshold_desc = f">= {LARGE_THRESHOLD}"

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
        array_setup=array_setup,
        index_calculation=index_calculation,
        country_lookup=country_lookup,
        storm_set_lookup=storm_set_lookup,
        sample_lookup=sample_lookup,
    )

    # Write script
    output_path.write_text(script_content)
    output_path.chmod(0o755)

    logger.info(f"Generated {output_path.name}:")
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
    parser.add_argument(
        "--no-check-missing",
        action="store_true",
        help="Disable checking for missing targets (generate all combinations)",
    )

    args = parser.parse_args()

    # Parse samples
    samples = [int(s.strip()) for s in args.samples.split(",")]
    logger.info(f"Samples to process: {samples}")
    logger.info(f"Check missing targets: {not args.no_check_missing}")
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
            category,
            categories[category],
            samples,
            storm_sets,
            output_path,
            check_missing=not args.no_check_missing,
        )

    logger.info("Script generation complete!")


if __name__ == "__main__":
    main()
