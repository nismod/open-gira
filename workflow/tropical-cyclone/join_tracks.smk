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
    resources:
        mem_mb = 48_000
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
