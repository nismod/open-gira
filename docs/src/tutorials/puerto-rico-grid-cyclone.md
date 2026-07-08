# Cyclone risk to Puerto Rico's grid

This tutorial follows the electricity/tropical-cyclone workflow: it estimates
how Hurricane Irma (2017) disrupted power supply in Puerto Rico. It assumes you
have finished [the Wales tutorial](wales-roads-flooding.md), or are otherwise
comfortable requesting result files from snakemake.

The shape of this workflow differs from flooding in one important way: it is
**event-based**. Instead of return-period hazard maps, you work with individual
storm tracks, estimate the wind field each one produces, fail grid components
above a wind-speed threshold, and re-allocate power across the degraded network
to see who loses supply.

> **On timings and sizes.** As in the first tutorial, the figures here are rough
> guides, and downloads dominate. Use `-n` dry runs to preview each step.
>
> *(maintainer to-do: add measured runtimes and disk usage, and confirm the
> storm and country identifiers below against the current data.)*

Throughout, we use one country and one storm:

- **Country:** `PRI` (Puerto Rico)
- **Storm:** `2017242N16333` (Irma, 2017), from the **IBTrACS** historic
  cyclone dataset

Later steps trigger the earlier steps they depend on, so you do not have to run
these strictly in order — but reading them in order shows the pipeline.

## Step 1 — Which storms hit which countries?

```bash
snakemake -n --cores 1 -- results/power/by_storm_set/IBTrACS/storms_by_country_impacted.json
```

This identifies, for the IBTrACS set, which storms are likely to have impacted
which countries. Drop the `-n` to run it.

## Step 2 — Storm tracks and wind fields

Preprocess the storm tracks into a common event-set format:

```bash
snakemake --cores 1 -- results/storm_tracks/IBTrACS/tracks.geoparquet
```

Estimate surface roughness over the grid region (from land-surface data), then
the maximum wind field each nearby storm produces, using a modified Holland wind
model downscaled to the surface:

```bash
snakemake --cores 1 -- results/power/by_country/PRI/storms/surface_roughness.tiff
snakemake --cores 1 -- results/power/by_country/PRI/storms/IBTrACS/max_wind_field.nc
```

## Step 3 — Expose the grid and re-allocate power

Split the grid edges on the wind grid, then, for our single storm, remove the
edges that exceed the configured wind-speed thresholds and re-allocate power from
plants to targets across the degraded network:

```bash
snakemake --cores 1 -- results/power/by_country/PRI/exposure/edges_split.geoparquet
snakemake --cores 1 -- results/power/by_country/PRI/exposure/IBTrACS/2017242N16333.nc
```

The wind-speed thresholds at which edges fail are set in `config/config.yaml`
under `transmission_windspeed_failure`; see the
[grid/cyclone usage page](../user-guide/usage/risk-analysis/grid-cyclone.md) for
the configuration.

## Step 4 — Aggregate the disruption

Run the exposure calculation for the storm across every country its track
crosses:

```bash
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/exposure_by_target.nc
```

## Step 5 — Visualise it

```bash
# Maximum wind speed over the impacted area
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/wind_field.png

# How supply degrades as the failure threshold varies
snakemake --cores 1 -- results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/outage_map/outage_map_by_threshold.gif
```

The result is a picture of Irma's wind field over Puerto Rico and an animation of
how modelled power supply degrades as the assumed failure threshold changes.

## Where to go next

- Swap `PRI` and the storm ID to study a different country or event.
- Replace `IBTrACS` (historic storms) with a synthetic storm set (STORM, IRIS,
  CHAZ) to explore return-period rather than single-event risk. The available
  storm sets and their scenario naming are in
  [Targets and wildcards](../reference/targets-and-wildcards.md).
- Read the [grid/cyclone usage page](../user-guide/usage/risk-analysis/grid-cyclone.md)
  for the full set of outputs and configuration.
