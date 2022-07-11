# Population data

Population data is obtained from https://www.worldpop.org/geodata/listing?id=79 and is used to allocate the appropriate GDP for target locations of interest as well as for statistical purposes.

Each country has a `{country}_ppp_2020_UNadj_constrained.tif` raster file which contains the geospatial population density. The `{country}` iso3 wildcard is retrieved from https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m. Internet connection is therefore required. The `download_population_all` will download the data for each country.