# Transport / flooding

The flooding risk analysis pipeline starts by creating an infrastructure
network (road or rail) as described previously. Please refer to this section to
configure the network creation.

## Description

The pipeline has the following core steps:

1. Download the raster files for the hazard dataset.
1. Crop the raster files to the extent of the network in question.
1. Split the network edges on the raster grid, so no edge resides in more than one pixel.
1. Find the pixel value / hazard intensity (e.g. flood depth) for every hazard raster, for each split edge.
1. For each asset type (e.g. unpaved road) where a damage curve is defined,
calculate the damage fraction for each split edge, for each raster. 
1. Using rehabilitation cost estimates, calculate the expected monetary loss due to the calculated damage fraction.
1. Concatenate the previously split edges, summing the expected monetary loss.
1. Given hazard rasters at various return periods, calculate the Expected Annual Damages (EAD).
1. Aggregate EAD to desired regional level.

## Configuration

The hazard component of the analysis is configurable in `config/config.yaml`:
- `hazard_datasets` contains hazard names pointing to files of hazard layers.
  These layers are currently flood rasters (inundation depths for a given return
  period). Add or amend an entry pointing to file containing the rasters you
  wish to consider.
- Ensure `hazard_types` contains an identical key referencing the hazard types.
  This is currently limited to `flood` only.
- Configure the damage curves:
    - Check and amend `direct_damages.asset_types` contains any assets you wish
      to calcuate direct damage costs for. Currently implemented assets are
      available in `src/open_gira/assets.py` as the classes inheriting from
      `Assets`.
    - Ensure `direct_damages.curves_dir` is set to a path containing damages
      functions per asset type, organised by hazard type. See
      `bundled_data/damage_curves` for an example.
    - These damage function files should be named `<asset_type>.csv`, e.g.
      `road_unpaved.csv`. Their format is exemplified by
      `bundled_data/damage_curves/flood/road_unpaved.csv`

## Outputs

### Rasterised network (per slice)

To generate a network, split by the hazard raster grid we might issue something like:
```bash
snakemake --cores all -- results/splits/<dataset_name>_filter-<network_type>/hazard-<hazard_name>/slice-<slice_number>.geoparquet
```

For the `egypt-latest` OSM dataset, `road` network, `aqueduct-river` hazard set and slice `0`, that would be:
```bash
snakemake --cores all -- results/splits/egypt-latest_filter-road/hazard-aqueduct-river/slice-0.geoparquet
```

### Expected Annual Damages (per slice)

To request an evaluation of Expected Annual Damages (EAD) as a function of
hazard Return Period (RP) for a given `slice`, we can request something like:

For example (with a `config.slice_count` of 9):
```bash
snakemake --cores all -- results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/slice-5.geoparquet
```

###  Expected Annual Damages (per admin region)

And to compute all the slices for a given domain and then aggregate to country level (admin level 0):
```bash
snakemake --cores all -- results/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet
```

For more possible outputs please refer to the detailed documentation and the
rules defined in `workflow/rules/`.