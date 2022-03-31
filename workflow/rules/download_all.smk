"""Download all files

Run all of the download rules, whether or not any futher rules in the pipeline require the data.
"""


rule download_all:
    input:
        out_fixed,
        out_events,
        out_population,
        out_GDP,
        out_powerplant,
        out_gridfinder,
        f"{config['output_dir']}/input/admin-boundaries/gadm36.gpkg",
        out_adminboundaries_codes,
