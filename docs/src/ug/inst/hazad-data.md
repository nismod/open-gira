# Tanzania flood hazard data

The last thing we need before we can get going is some data on the hazards that threaten the infrastructure.
The hazard information is presented as a series of image files.
Each image file shows the expected flood coverage within a given return 
period.[^return period]
The images are also split by whether they include subsidence in the modelling, the climate model they use,
and when their data or predictions come from.

Although the .tif files are image files, they do not have normal image information (so they will look completely black
if opened in an image editor). 
Instead of several colour/alpha channels of data, these files just have a single channel of data with (usually) very
low values (specifying water height at each pixel).
To see the flood information in an image, open it in QGIS.
In the screenshot below, we chose the last file in the `./data/aqueduct/` directory,
`./data/aqueduct/inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000.tif`.
Alongside the image (.tif) files, there are two .csv files that include metadata for the hazard data.
These .csv files (one for coastal and one for river flood hazards) describe what each image file represents.

![QGIS screenshot showing flood height data.](../../img/QGIS-hazard.png)

You can acquire the .zip file with the images and metadata from 
[the Zenodo repository](https://zenodo.org/record/5887564).
All these files should be unpacked and placed in the `./data/aqueduct` directory, so that the path to 
aqueduct_river.csv is `./data/aqueduct/aqueduct_river.csv`, and all the .csv and .tif files are in the same directory.

[^return period]: This means the most dramatic flooding we would expect within a given window of time, e.g. 10 years or
    100 years.
