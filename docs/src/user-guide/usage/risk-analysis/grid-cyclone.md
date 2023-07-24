# Electricity grids - tropical cyclones

This analysis intersects electricity grids with tropical cyclone wind speed
risk. Network creation is described previously. The hazards are event-based (as
opposed to return period map based).

## Description

1. Download storm track data.
1. Process storm track data into common geoparquet format, ensure meterological variables are consistent.
1. Define a wind speed grid, with an extent equal to the network bounding box, plus a buffer.
1. Given a storm track, estimate the resulting wind speed field as a function of time using a modified Holland model.
1. Downscale the wind speeds from free-atmosphere to surface level, using a power law. The downscaling exponent is a function of the surface roughnesses, derived from a land surface categorisation.
1. Take the maximum wind speed across time for each pixel in the wind grid.
1. For a certain defined threshold or set of thresholds, we then fail electricity grid edges (transmission and distribution lines) that experience
wind in excess of the threshold.
1. Allocate power across the newly degraded grid and see how it differs to the nominal allocation.
1. Report on change in power supply and numbers of people affected.

## Configuration

- Review and amend `config/config.yaml`:
    - `storm_sets` should contain a storm set name, pointing to a JSON file.
      This JSON file should contain an empty list to process all storms for
      this storm set, or a list of string storm IDs if only a subset is
      required.
    - Specify `transmission_windspeed_failure` as a list of wind speeds in
      meters per second at which to fail edges.

See comments in the `config/config.yaml` file for other less crucial
configuration options.

## Outputs

Below you will find example commands necessary to create outputs for the case of
__PRI__ (Puerto Rico) and the storm __2017242N16333__ (Irma, 2017) from the
__IBTrACS__ historic cyclone dataset. Later steps will trigger earlier,
dependent steps if necessary.

### Networks at risk
- Identify which storms are likely to impact which countries
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/storms_by_country_impacted.json`

### Maximum wind speed fields
- Download storm track data (historic or synthetic)
- Preprocess data into a (mostly) common event-set format.
    - `snakemake --cores 1 -- results/storm_tracks/IBTrACS/tracks.geoparquet`
- Download land surface categorical data
- Estimate surface roughness for the grid region
    - `snakemake --cores 1 -- results/power/by_country/PRI/storms/surface_roughness.tiff`
- Any storm tracks within some proximity of the grid in question are processed into maximum wind speed footprints for the electricity grid region. This is by means of a modified Holland wind field model. Downscale the winds to the surface by a power law, using the surface roughness data.
    - `snakemake --cores 1 -- results/power/by_country/PRI/storms/IBTrACS/max_wind_field.nc`

## Network exposure
- Rasterise the electricity grid (split line segments such that no segment extends beyond a wind cell boundary).
    - `snakemake --cores 1 -- results/power/by_country/PRI/exposure/edges_split.geoparquet`
- For a given storm and country, remove edges from the grid which exceed a given wind speed damage threshold. Reallocate power from powerplants to targets, by GDP, on the degraded network. Store the ratio between the original power supplied and the degraded power supplied to each target.
    - `snakemake --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/2017242N16333.nc`
- Perform the above calculations for a single storm in a cyclone dataset and all the countries its track intersects.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/exposure_by_target.nc`
- Perform the above calculations for an entire cyclone dataset and all the countries its tracks intersect.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/exposure.txt`

## Further analysis

TODO: Add additional rules to facilitate comparison across countries and cyclone datasets.

TODO: Document these additional rules.

- Map how electricity supply is impacted by a given storm for a variety of different damage thresholds.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/outage_map/outage_map_by_threshold.gif`
- Map the maximum wind speed experienced across the area impacted by a storm.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/wind_field.png`

There also exist other plotting and mapping steps to visualise intermediate and final outputs. Refer to `workflow/rules/analysis` for a description of these.