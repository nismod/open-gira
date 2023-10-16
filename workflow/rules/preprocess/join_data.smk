rule join_data:
    """
    Take split .geoparquet files and output a single, unified .geoparquet file
    for each dataset-hazard combination
    """
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
snakemake --cores all results/tanzania-mini_filter-road_hazard-aqueduct-river.geoparquet
"""


def all_storm_tracks_by_sample(wildcards) -> list[str]:
    """
    Return a list of every per-sample tracks file for a given STORM_SET.
    """
    dataset_name = wildcards.STORM_SET.split("-")[0]
    return expand(
        "{OUTPUT_DIR}/storm_tracks/{STORM_SET}/{SAMPLE}/tracks.geoparquet",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        STORM_SET=wildcards.STORM_SET,
        SAMPLE=range(0, SAMPLES_PER_TRACKSET[dataset_name])
    )


rule concat_storm_tracks:
    """
    To concat all samples of tracks together, useful when we want to e.g. plot track atlas.
    """
    input:
        by_sample=all_storm_tracks_by_sample
    output:
        tracks_from_all_samples="{OUTPUT_DIR}/storm_tracks/{STORM_SET}/tracks.geoparquet",
    run:
        import geopandas as gpd
        import pandas as pd

        data = [gpd.read_parquet(path) for path in input.by_sample]
        concat = gpd.GeoDataFrame(pd.concat(data))
        concat.to_parquet(output.tracks_from_all_samples)

"""
Test with:
snakemake -c1 -- results/storm_tracks/STORM-constant/tracks.geoparquet
"""