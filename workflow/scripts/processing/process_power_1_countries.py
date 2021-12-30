"""
This file downloads the plants and targets data to csv files. Can be made to run in parallel.
"""

from process_power_functions import write_plants_targets
import sys

# Note: snakemake will find if file exists already
write_plants_targets(sys.argv[1], sys.argv[2], sys.argv[3])
