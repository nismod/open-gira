# Administrative boundaries

The administrative boundaries are downloaded from https://gadm.org/data.html. This data is needed to determine the appropriate GDP value for a location of interest, cross-checking country data and also for statistical purposes. There are three rules which may be executed.

The `download_gadm` rule will download the 'gadm36_gpkg.zip' file which contains the administrative regions globally.

The `download_gadm_by_country` rule will download the individually selected country data of 'gadm36_gpkg.zip' file (see above).

Lastly, the `download_gadm_levels` rule will download the equivalent of 'gadm36_gpkg.zip' but as five separate layers, where each is a level of subdivision/aggregation (ranging from smallest administrative region to country). See https://gadm.org/download_world.html for further details.

