# open-gira configuration directory

## README.md
This file.

## config.yaml
YAML configuration file. This file contains parameters that are frequently
modified prior to an analysis. Refer to it for more information.

## land_use_to_surface_roughness.csv
A lookup table from a land use classification scheme to a typical surface
roughness length. This is used in downscaling wind speeds in tropical-cyclone
analyses.

## arc_cluster
`snakemake` profile for running open-gira on the University of Oxford's ARC
cluster.

## damage_curves
CSV files, one per `asset_type`, describing the relationship between hazard
intensity and damage fraction. These are used for direct damage estimation.

## hazard_resource_locations
Fairly self-explanatory, contains files containing lists of hazard file URLs to
acquire for analyses.

## osm_filters
Text files for use by osmium to [filter features based on
tag](https://docs.osmcode.org/osmium/latest/osmium-tags-filter.html). We used
these to construct infrastructure networks from OpenStreetMap data.

## rehab_costs
CSV files containing rehabilitation cost estimates in USD. One file per sector,
each row a unique `asset_type`.

## storm_sets
JSON files, each containing a single list. The list may contain zero or many
unique string storm IDs. These can be used to subset a parent `STORM_SET`, e.g.
IBTrACS may be subset to contain just storms that hit Mexico since 2000, if the
relevant storm IDs are provided. Where the list provided is empty, this is
interpreted as "process every storm in the parent set", i.e. all of IBTrACS.
