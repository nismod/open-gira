# Download data

The analysis requires the following data
- [Administrative boundaries](adminboundaries.md)
- [GDP data](GDP.md)
- [Gridfinder data](gridfinder.md)
- [Population data](population.md)
- [Powerplant data](powerplants.md)
- [STORM data](stormtracks.md)

This data is not interlinked (order irrelevant) and forms the following simple workflow diagram for the `download_all` rule.

![`download_all` rule workflow.](../../img/dag_download_all.png)

This data can be downloaded in one by entering `snakemake -s workflow/Snakefile download_all -c{core_numbers}` in the linux command line when in the root open-gira folder. Note that each individual bullet link explains how to download individual data selections, however, it is recommended to use the aforementioned command as snakemake will automatically detect already downloaded files in one go.
