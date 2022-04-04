# Configuration

The configuration file for the run is found in `./config/config.yaml`.

Your configuration file should look like the below:

```yaml
# Output directory
output_dir: 'results'

# Named list of OpenStreetMap datasets
infrastructure_datasets:
  tanzania-latest: 'https://download.geofabrik.de/africa/tanzania-latest.osm.pbf'
  wales-latest: 'https://download.geofabrik.de/europe/great-britain/wales-latest.osm.pbf'

# Named list of .txt files that specify hazard raster file locations
hazard_datasets:
  aqueduct-coast: 'https://raw.githubusercontent.com/mjaquiery/aqueduct/main/tiffs.txt'
  aqueduct-river: 'https://raw.githubusercontent.com/mjaquiery/aqueduct/main/rivers.txt'

# Number of slices to cut dataset into -- must be a square number
slice_count: 36
# Edge attributes to preserve during network/hazard intersection.
keep_tags: 'highway'
# Filters definition -- base filename is used as the filter name
osmium_tags_filters_file: "config/highway-core.txt"

# Options for exposure tif file generation
exposure_tifs:
  # When generating exposure summaries, filter by cells with at >= exposure_threshold m of flooding
  exposure_threshold: 0.5
  # Used when scaling raster files, values < 1 will downsample
  scaling_factor: 0.1
  # Used to determine how resampling occurs. Valid values are in rasterio.enums.Resampling.
  resampling_mode: 'bilinear'
  # Options for plotting
  plot:
    # Keyword arguments for raster plotting
    raster:
      # Colour mapping to use for raster data
      cmap: 'Reds'
```

You should also have a `./config/highway-core.txt` file that looks like:
```text
motorway,motorway_link,trunk,trunk_link,primary,primary_link,secondary,secondary_link
```

This configuration file holds everything we need to run the workflow.
The directory `<output_dir>` (hereafter `./results/`) will be created when the workflow runs.
The infrastructure datasets will be downloaded, along with all the hazard raster 
files listed in the hazard datasets.
The infrastructure datasets will be filtered and sliced, and then combined with
the hazard data to produce four output files, one for each of the 
infrastructure-hazard dataset combinations.

The infrastructure datasets listed in this config are real, but the hazard datasets
only contain a tiny fraction of the raster files available. 
This means the workflow is relatively quick to run compared to 
when we have gigabytes of hazard raster files.
