rule parse_ibtracs:
    input:
        ibtracs_csv = "results/input/IBTrACS/raw/v4.csv"
    output:
        ibtracs_parquet = "results/input/IBTrACS/processed/v4.parquet"
    script:
        "../../scripts/process/parse_IBTrACS.py"

"""
To test:
snakemake -c1 results/input/IBTrACS/processed/v4.parquet
"""


rule slice_ibtracs:
    input:
        ibtracs=rules.parse_ibtracs.output.ibtracs_parquet,
        global_boxes=rules.world_splitter.output.global_boxes,
    output:
        sliced_tracks="{OUTPUT_DIR}/power/slice/{BOX}/storms/IBTrACS.geoparquet",
    run:
        import pandas as pd
        import geopandas as gpd

        boxes = gpd.read_parquet(input.global_boxes).set_index("box_id")
        box_geom = boxes.loc[f"box_{wildcards.BOX}", "geometry"]

        df = pd.read_parquet(input.ibtracs)
        tracks = gpd.GeoDataFrame(
            data=df,
            geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
        )
        tracks[tracks.intersects(box_geom)].to_parquet(output.sliced_tracks)


"""
To test:
snakemake -c1 results/power/slice/1030/storms/IBTrACS.geoparquet
"""
