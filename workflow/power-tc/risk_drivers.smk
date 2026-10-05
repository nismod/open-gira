"""
Collate tables of outage event data, along with storm and network characteristics.
"""


def pop_affected_by_country_paths(wildcards) -> list[str]:
    import json
    country_set_path = checkpoints.countries_intersecting_storm_set.get(
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        STORM_SET=wildcards.STORM_SET,
    ).output.country_set
    with open(country_set_path, "r") as fp:
        countries = json.load(fp)
    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/disruption/{STORM_SET}/pop_affected_by_event.pq",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=countries,
        STORM_SET=wildcards.STORM_SET,
    )

def tracks_by_country_paths(wildcards) -> list[str]:
    import json
    country_set_path = checkpoints.countries_intersecting_storm_set.get(
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        STORM_SET=wildcards.STORM_SET,
    ).output.country_set
    with open(country_set_path, "r") as fp:
        countries = json.load(fp)
    return expand(
        "{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{STORM_SET}/0/tracks.geoparquet",
        OUTPUT_DIR=wildcards.OUTPUT_DIR,
        COUNTRY_ISO_A3=countries,
        STORM_SET=wildcards.STORM_SET,
    )

rule risk_drivers_by_country_by_storm_set:
    """
    Test with:
    snakemake -c1 results/power/by_storm_set/IBTrACS/disruption/drivers.pq
    """
    input:
        all_disruption = pop_affected_by_country_paths,
        all_tracks = tracks_by_country_paths,
    output:
        drivers = "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/disruption/drivers.pq"
    run:
        import geopandas as gpd
        import pandas as pd

        countries = [p.split("/")[-4] for p in input.all_disruption]
        data = []
        for iso_a3 in countries:
            disruption = pd.read_parquet(f"{wildcards.OUTPUT_DIR}/power/by_country/{iso_a3}/disruption/{wildcards.STORM_SET}/pop_affected_by_event.pq")
            disruption.columns = "pop_" + disruption.columns
            track = gpd.read_parquet(f"{wildcards.OUTPUT_DIR}/power/by_country/{iso_a3}/storms/{wildcards.STORM_SET}/0/tracks.geoparquet")
            if disruption.empty or track.empty:
                continue
            else:
                extremes = track.groupby("track_id").agg(
                    min_pressure_hpa=("min_pressure_hpa", "min"),
                    max_wind_speed_ms=("max_wind_speed_ms", "max"),
                )
                join = disruption.join(extremes).sort_index().reset_index()
                join["iso_a3"] = iso_a3
                join["storm_set"] = wildcards.STORM_SET
            data.append(join)
        df = pd.concat(data)
        # Put pop_ columns at the end
        df = df[sorted(df.columns, key=lambda col: col.startswith("pop_"))]
        df.to_parquet(output.drivers)

