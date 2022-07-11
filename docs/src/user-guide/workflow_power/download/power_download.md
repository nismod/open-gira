# Download data for power analysis

The power analysis requires the following data
- [Administrative boundaries](power_download_adminboundaries.md)
- [GDP data](power_download_GDP.md)
- [Gridfinder data](power_download_gridfinder.md)
- [Population data](power_download_population.md)
- [Powerplant data](power_download_powerplants.md)
- [STORM data](power_download_stormtracks.md)


This data is not interlinked (order irrelevant) and forms the following simple workflow diagram for the `download_all` rule.

![`download_all` rule workflow.](../power_img/dag_download_all.png)

This data can be downloaded in one by entering `snakemake -s workflow/Snakefile download_all -c{core_numbers}` in the linux command line when in the root open-gira folder. Note that each individual bullet link explains how to download individual data selections, however, it is recommended to use the aforementioned command as snakemake will automatically detect already downloaded files in one go.


