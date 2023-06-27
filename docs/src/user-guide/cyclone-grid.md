# Electricity grids - tropical cyclones

This workflow estimates how tropical cyclones impact upon electricity grids. The spatial domain can range from country-specific to global.

The major steps are as follows. Some of the steps are given with example commands necessary to invoke them for the case of __PRI__ (Puerto Rico) and the storm __2017242N16333__ (Irma, 2017) from the __IBTrACS__ historic cyclone dataset. Later steps will trigger earlier, dependent steps if necessary.

### Determine which networks are at risk
- Identify which storms are likely to impact which countries
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/storms_by_country_impacted.json`

### Build electricity networks
- Download gridfinder electricity grid data (line edges, 'target' consuming polygons)
- Create global 'targets' file, where each electricity consuming target is annotated with the population and GDP of the enclosed area.
    - `snakemake --cores 1 -- results/power/targets.geoparquet`
- Download WRI powerplant data
- Download GDP per capita data
- Download population data
- Construct an electricity grid for each country to be studied where nodes are power generators or target consumers. Power consumption is proportional to the GDP enclosed by each target polygon.
    - `snakemake --cores 1 -- results/power/by_country/PRI/network/edges.geoparuqet`

### Generate maximum wind speed fields
- Download storm track data (historic or synthetic)
- Preprocess data into a (mostly) common event-set format.
    - `snakemake --cores 1 -- results/storm_tracks/IBTrACS/tracks.geoparquet`
- Download land surface categorical data
- Estimate surface roughness for the grid region
    - `snakemake --cores 1 -- results/power/by_country/PRI/storms/surface_roughness.tiff`
- Any storm tracks within some proximity of the grid in question are processed into maximum wind speed footprints for the electricity grid region. This is by means of a modified Holland wind field model. Downscale the winds to the surface by a power law, using the surface roughness data.
    - `snakemake --cores 1 -- results/power/by_country/PRI/storms/IBTrACS/max_wind_field.nc`

### Expose networks to maximum wind speeds
- Rasterise the electricity grid (split line segments such that no segment extends beyond a wind cell boundary).
    - `snakemake --cores 1 -- results/power/by_country/PRI/exposure/edges_split.geoparquet`
- For a given storm and country, remove edges from the grid which exceed a given wind speed damage threshold. Reallocate power from powerplants to targets, by GDP, on the degraded network. Store the ratio between the original power supplied and the degraded power supplied to each target.
    - `snakemake --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/2017242N16333.nc`
- Perform the above calculations for a single storm in a cyclone dataset and all the countries its track intersects.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/exposure_by_target.nc`
- Perform the above calculations for an entire cyclone dataset and all the countries its tracks intersect.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/exposure.txt`

### Analyse exposure results

TODO: Add additional rules to facilitate comparison across countries and cyclone datasets.

TODO: Document these additional rules.

- Map how electricity supply is impacted by a given storm for a variety of different damage thresholds.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/outage_map/outage_map_by_threshold.gif`
- Map the maximum wind speed experienced across the area impacted by a storm.
    - `snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/wind_field.png`

There also exist other plotting and mapping steps to visualise intermediate and final outputs. Refer to `workflow/rules/analysis` for a description of these.
