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
