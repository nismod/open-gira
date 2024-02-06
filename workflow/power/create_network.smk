"""
Network creation routines
"""

import pandas as pd

from open_gira.io import cached_json_file_read
from open_gira.curves import logistic_min


rule gridfinder_to_geoparquet:
    """
    Store linestrings as geoparquet for faster accessing.
    Store a representative point from each linestring for faster spatial joins.
    """
    input:
        geopackage = rules.download_gridfinder.output.electricity_grid_global,
    resources:
        mem_mb = 4096
    output:
        linestring = "{OUTPUT_DIR}/power/gridfinder.geoparquet",
        rep_point = "{OUTPUT_DIR}/power/gridfinder_rep_point.geoparquet",
        plot = "{OUTPUT_DIR}/power/gridfinder.png",
    run:
        import geopandas as gpd
        import pandas as pd
        import datashader as ds
        import spatialpandas
        import spatialpandas.io

        grid = gpd.read_file(input.geopackage)

        # write out linestrings as geoparquet
        grid.to_parquet(output.linestring)

        # cast to spatialpandas for plotting
        grid_sp = spatialpandas.GeoDataFrame(grid)

        # make an integer source category column
        cat = pd.get_dummies(grid_sp.source)
        cat.gridfinder *= 2
        # openstreetmap = 1, gridfinder = 2
        grid_sp["source_id"] = cat.sum(axis=1).astype(int)

        # plot the gridfinder network
        cvs = ds.Canvas(plot_width=1500, plot_height=670)
        agg = cvs.line(grid_sp, geometry='geometry', agg=ds.mean("source_id"))
        img = ds.transfer_functions.shade(agg)
        ds.utils.export_image(img=img, filename=output.plot.split(".")[0], fmt=".png", background="black")

        # pick a point somewhere on each linestring, replace geometry with this
        grid["geometry"] = grid.geometry.representative_point()

        # write out representative points
        spatialpandas.io.to_parquet(spatialpandas.GeoDataFrame(grid), output.rep_point)

"""
To test:
snakemake --cores 1 results/input/gridfinder/gridfinder.geoparquet
"""


rule subset_gridfinder:
    """
    Subset the gridfinder dataset to a country boundary. Can be quite a heavy
    operation depending on the size of the country.

    N.B. Should take around an hour for the USA/gridfinder case. Could be
    reduced by using dask workers from spatialpandas, but probably not worth
    the complexity?
    """
    input:
        gridfinder=rules.gridfinder_to_geoparquet.output.linestring,
        gridfinder_rep_point=rules.gridfinder_to_geoparquet.output.rep_point,
        admin_bounds="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet"
    output:
        gridfinder="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/gridfinder.geoparquet",
    resources:
        mem_mb=8192
    run:
        import geopandas as gpd
        import spatialpandas
        import spatialpandas.io

        import logging

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        os.makedirs(os.path.dirname(output.gridfinder), exist_ok=True)

        # read admin bounds for country in question
        countries: gpd.GeoDataFrame = gpd.read_parquet(input.admin_bounds)
        country_gp: gpd.GeoDataFrame = countries[countries.GID_0 == wildcards.COUNTRY_ISO_A3]

        grid: gpd.GeoDataFrame = gpd.read_parquet(input.gridfinder)

        if country_gp.geometry.area.sum() > 10:  # square degrees, not a fair measure at high latitudes

            logging.info(f"Using spatialpandas point-in-polygon to subset gridfinder for {wildcards.COUNTRY_ISO_A3}")
            # create a spatialpandas GeoDataFrame for the country
            # it tends to be faster than geopandas for sjoins, about 3x in Mexico test case
            country_sp: spatialpandas.GeoDataFrame = spatialpandas.GeoDataFrame(country_gp)

            # read in representative points of linestrings as spatialpandas geodataframe
            # spatialpandas sjoin can only do point-in-polygon, not linestring-in-polygon
            grid_rep_point: spatialpandas.GeoDataFrame = spatialpandas.io.read_parquet(input.gridfinder_rep_point)

            logging.info(f"Spatially joining gridfinder with {wildcards.COUNTRY_ISO_A3}")
            grid_subset: spatialpandas.GeoDataFrame = spatialpandas.sjoin(grid_rep_point, country_sp, how="inner")

            logging.info(f"Writing gridfinder to {output.gridfinder}")
            # use the joined index to select from the geopandas linestring geodataframe
            grid.loc[grid_subset.index, :].to_parquet(output.gridfinder)

        else:

            logging.info(f"Spatially joining gridfinder with {wildcards.COUNTRY_ISO_A3}")
            grid_subset: gpd.GeoDataFrame = grid.sjoin(country_gp, how="inner")

            logging.info(f"Writing gridfinder to {output.gridfinder}")
            grid_subset[grid.columns].to_parquet(output.gridfinder)

"""
Test with:
snakemake -c1 results/power/by_country/HTI/network/gridfinder.geoparquet
"""


rule subset_targets:
    """
    Subset the targets dataset to a country boundary
    """
    input:
        targets=rules.annotate_targets.output.targets,
        admin_bounds="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet"
    resources:
        mem_mb=8192
    output:
        targets="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
    run:
        import geopandas as gpd
        import spatialpandas

        import logging

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        os.makedirs(os.path.dirname(output.targets), exist_ok=True)

        countries: gpd.GeoDataFrame = gpd.read_parquet(input.admin_bounds)
        country_gp: gpd.GeoDataFrame = countries[countries.GID_0 == wildcards.COUNTRY_ISO_A3]

        targets: gpd.GeoDataFrame = gpd.read_parquet(input.targets)

        if country_gp.geometry.area.sum() > 10:  # square degrees, not a fair measure at high latitudes

            logging.info(f"Using spatialpandas point-in-polygon to subset targets for {wildcards.COUNTRY_ISO_A3}")
            # create a spatialpandas GeoDataFrame for the country
            # it tends to be faster than geopandas for sjoins, about 3x in Mexico test case
            country_sp: spatialpandas.GeoDataFrame = spatialpandas.GeoDataFrame(country_gp)

            # spatialpandas sjoin can only do point-in-polygon, not linestring-in-polygon
            # so make some representative points
            targets_rep_point = targets.copy()
            targets_rep_point.geometry = targets_rep_point.geometry.representative_point()
            targets_rep_point = spatialpandas.GeoDataFrame(targets_rep_point)

            logging.info(f"Spatially joining targets with {wildcards.COUNTRY_ISO_A3}")
            targets_subset: spatialpandas.GeoDataFrame = spatialpandas.sjoin(targets_rep_point, country_sp, how="inner")

            logging.info(f"Writing targets to {output.targets}")
            # use the joined index to select from the geopandas linestring geodataframe
            targets.loc[targets_subset.index, :].to_parquet(output.targets)

        else:

            logging.info(f"Spatially joining targets with {wildcards.COUNTRY_ISO_A3}")
            targets_subset: gpd.GeoDataFrame = targets.sjoin(country_gp, how="inner")

            logging.info(f"Writing targets to {output.targets}")
            targets_subset[targets.columns].to_parquet(output.targets)

"""
Test with:
snakemake -c1 results/power/by_country/HTI/network/targets.geoparquet
"""


rule subset_powerplants:
    """
    Subset the powerplants dataset to a country boundary
    """
    input:
        admin_bounds="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
        powerplants="{OUTPUT_DIR}/power/powerplants.geoparquet",
    output:
        powerplants="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/powerplants.geoparquet",
    run:
        import geopandas as gpd

        import logging

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        os.makedirs(os.path.dirname(output.powerplants), exist_ok=True)

        countries: gpd.GeoDataFrame = gpd.read_parquet(input.admin_bounds)
        country_gp: gpd.GeoDataFrame = countries[countries.GID_0 == wildcards.COUNTRY_ISO_A3]

        powerplants: gpd.GeoDataFrame = gpd.read_parquet(input.powerplants)

        logging.info(f"Spatially joining powerplants with {wildcards.COUNTRY_ISO_A3}")
        powerplants_subset: gpd.GeoDataFrame = powerplants.sjoin(country_gp, how="inner")

        logging.info(f"Writing targets to {output.powerplants}")
        powerplants_subset[powerplants.columns].to_parquet(output.powerplants)

"""
Test with:
snakemake -c1 results/power/by_country/HTI/network/powerplants.geoparquet
"""


def threads_for_country(wildcards) -> int:
    """
    Estimate a sensible number of threads to use for a given country. We rank
    countries by their number of targets, and then apply a sigmoid (logistic
    minimum) function to the ranking.

    N.B. Rules that employ this function must also include
    `country_target_count_path` as an input.

    Args:
        wildcards: Must include COUNTRY_ISO_A3 to do the country lookup.

    Returns:
        Thread allocation
    """

    ranking_file = checkpoints.rank_countries_by_target_count.get(**wildcards).output.lookup_table
    ranked = pd.read_csv(ranking_file)

    ranked["threads"] = logistic_min(
        ranked.index,  # input to transform
        0.85 * workflow.cores,  # maximum (roughly)
        -2,  # minimum (roughly)
        0.02,  # steepness of sigmoid
        0.8 * len(ranked.index)  # location of sigmoid centre on input axis
    ).astype(int)

    try:
        n_threads, = ranked.loc[ranked.iso_a3 == wildcards.COUNTRY_ISO_A3, "threads"]
    except ValueError:
        # likely country is missing from table
        n_threads = 1

    # ensure it's at least one
    return max([1, n_threads])


def country_target_count_path(wildcards) -> str:
    """
    We depend on the file (path) returned by this function to allocate
    resources (number of CPUs and therefore memory) to certain rules. Those
    rules use `threads_for_country` as a function returning a value for the
    `threads` parameter. However those rules must also include this function as
    an `input` function, to ensure the CSV data file is available for
    `threads_for_country` to execute. Having the checkpoint lookup within
    `threads_for_country` is not, on its own, sufficient.
    """
    return checkpoints.rank_countries_by_target_count.get(**wildcards).output.lookup_table


rule create_power_network:
    """
    Combine power plant, consumer and transmission data for given area
    """
    input:
        # `threads_for_country` will fail unless this CSV is present when resources are set
        country_target_count=country_target_count_path,
        plants="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/powerplants.geoparquet",
        targets="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/targets.geoparquet",
        gridfinder="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/gridfinder.geoparquet",
    threads: threads_for_country
    output:
        edges="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        nodes="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/nodes.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json",
    script:
        "./create_electricity_network.py"

"""
Test with:
snakemake -c1 results/power/by_country/HTI/edges.geoparquet
"""


def networks_affected_by_storm_set(wildcards):
    """
    Given STORM_SET as a wildcard, lookup the countries that a storm set
    affects and return paths to their network edge files.
    """

    json_file = checkpoints.countries_intersecting_storm_set.get(**wildcards).output.country_set
    country_set = cached_json_file_read(json_file)

    return expand(
        "results/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        COUNTRY_ISO_A3=country_set,  # list of str
        STORM_SET=wildcards.STORM_SET  # str
    )


rule create_networks_affected_by_storm_set:
    """
    A target rule to create all networks potentially impacted by a storm set.
    """
    input:
        networks = networks_affected_by_storm_set
    output:
        completion_flag = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/networks.txt"
    shell:
        """
        # one output file per line
        echo {input.networks} | tr ' ' '\n' > {output.completion_flag}
        """

"""
Test with:
snakemake -c1 -- results/power/by_storm_set/IBTrACS/networks.txt
"""


rule map_network_components:
    """
    Produce a datashader plot of edges within network, coloured by component ID.

    Useful for spotting electricity 'islands'.
    """
    input:
        edges="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
    output:
        plot="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.png",
    run:
        import geopandas as gpd
        import datashader as ds
        import datashader.transfer_functions as tf
        import matplotlib.cm
        import spatialpandas

        edges = gpd.read_parquet(input.edges)

        minx, miny, maxx, maxy = edges.total_bounds
        x = maxx - minx
        y = maxy - miny
        aspect = x / y
        width = 800
        cvs = ds.Canvas(plot_width=width, plot_height=int(width / aspect))

        agg = cvs.line(spatialpandas.GeoDataFrame(edges), geometry='geometry', agg=ds.mean("component_id"))
        img = tf.shade(agg, cmap=matplotlib.cm.Set2_r)

        ds.utils.export_image(img=img, filename=output.plot.split(".")[0], fmt=".png", background="black")

"""
Test with:
snakemake -c1 results/power/by_country/HTI/edges.png
"""
