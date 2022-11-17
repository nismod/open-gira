"""
File conversion rules
"""


rule geopackage_to_geoparquet:
    input:
        geopackage="{BASENAME}.gpkg",
    output:
        geoparquet="{BASENAME}.geoparquet",
    run:
        import geopandas as gpd
        gpd.read_file(input.geopackage).to_parquet(output.geoparquet)

"""
To test:
snakemake --cores 1 results/input/gridfinder/grid.geoparquet
"""
