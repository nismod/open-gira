# Configuration

The configuration file for the run is found in `./config/config.yaml`.
If everything has been set up according to this guide so far, this config file should not need changing.
If the infrastructure or hazard data are stored somewhere other than `./data` and `./data/aqueduct` then
the relevant paths should be specified in `data_dir` and `hazard_data_dir`.
The `output_dir` will be created when the workflow is run.
As with the other paths, this will be interpreted with respect to the open-gira directory, so the default
`results` will create `./results` to house the output.
Make sure dataset is `tanzania-latest`. 
Everything else should stay as it is.

## Setup

The final piece we need to have in place prior to running is a file that explains to open-gira how to slice up
the data files we gave it.
This file will be `./data/tanzania-latest-extracts.geojson` and we will generate it automatically by using
the `prepare-extracts.py` script.

First, however, we need to create a `./tanzania-latest.json` file that will be sliced into extracts by the script.
To do this, we create a new file, `./tanzania-latest.json` which looks like this:

```json
{
    "directory": "./data",
    "extracts": [
        {
            "bbox": [
                29.24395,
              -11.775945,
              40.69487,
              -0.974988
            ],
            "output": "tanzania-latest.osm.pbf"
        }
	]
}
```

Those values in the `bbox` property are the start and end latitude and longitude of the area we want to slice up.
We can extract the values for the whole `tazania-latest.osm.pbf` file by looking in the header information of the file.
Run:

```shell
osmium fileinfo data/tazania-latest.osm.pbf
```

You should see an output that contains
```text
Header:
  Bounding boxes:
    (29.24395,-11.775945,40.69487,-0.974988)
```

Once we've created that `./tanzania-latest.json` file, we can create the extracts definition file automatically 
by running:

```
python prepare-extracts.py tanzania-latest.json 6
```

The script creates an n-by-n grid where n is the second command line argument, 6 in the example above.
Once the command is run we should see it has created a file `./data/tanzania-latest-extracts.geojson` with 36 (i.e. 6^2) 
slices defined in it (0 through 35).