# Take split .geoparquet files and output a single, unified .geoparquet file
# for each dataset-hazard combination


rule join_exposure:
    input:
        lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "splits",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "slice-{i}.geoparquet",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
    script:
        "../../scripts/join_data.py"


"""
Test with:
snakemake --cores all results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet
"""


rule join_direct_damage_fraction:
    input:
        slices = lambda wildcards: expand(
            os.path.join(
                "{OUTPUT_DIR}",
                "direct_damages",
                "{DATASET}_{FILTER_SLUG}",
                "{HAZARD_SLUG}",
                "fraction",
                "slice-{i}.geoparquet",
            ),
            **wildcards,
            i=range(config["slice_count"])
        ),
    output:
        joined = "{OUTPUT_DIR}/direct_damages/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/damage_fraction.geoparquet",
    run:
        import os

        import pandas as pd

        parent_dir = os.path.dirname(input.slices[0])
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

        unified = gpd.GeoDataFrame(pd.concat([gpd.read_parquet(path) for path in input.slices]))
        unified.to_parquet(output.joined)


"""
Test with:
snakemake --cores all results/egypt-latest_filter-road_hazard-aqueduct-river_damage-fraction.parquet
"""
