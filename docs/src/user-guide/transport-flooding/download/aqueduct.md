# Hazard raster data download

The hazard raster files are not listed directly in the 
[config file](../configuration.md) because they tend to be numerous.
Instead, the config file points to a .txt file that lists where the
raster files can be obtained.
The name of the config list item specifies the hazard that the files represent.

First, the .txt file is downloaded (if it's on the internet) or copied.
Then, for each line in the text file, the raster file is downloaded (if it's on
the internet) or copied into `./results/input/hazard-<hazard>/raw/`.

For our first config `hazard_datasets` entry, `aqueduct-coast`, we specified the file
`https://raw.githubusercontent.com/mjaquiery/aqueduct/main/tiffs.txt`.
The contents of this file are:

```text
http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inuncoast_historical_nosub_hist_rp0001_5.tif
http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inuncoast_historical_nosub_hist_rp0002_0.tif
http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inuncoast_historical_nosub_hist_rp0005_0.tif
http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inuncoast_historical_nosub_hist_rp0010_0.tif
```

Each line gives the location of a .tif file that contains hazard data.
Each of these files will be downloaded and saved to `./results/input/hazard-aqueduct-coast/raw/`.