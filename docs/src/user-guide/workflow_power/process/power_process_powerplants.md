# Process power plants

In this file, the [power plant input data](../download/power_download_powerplants.md) is process to contain the columns
"source_id", "name", "type", "capacity_mw", "estimated_generation_gwh_2017", "primary_fuel", "box_id", "geometry" and
then split into a .csv file for each box_id. Note that if no power plants exist in the box_id, then an empty .csv file is
produced. This is to allow for consistent snakemake workflows as it can not be known in advance.
