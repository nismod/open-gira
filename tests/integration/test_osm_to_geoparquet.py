from . import runner


def test_osm_to_geoparquet():
    runner.run_snakemake_test(
        "osm_to_geoparquet",
        (
            "results/geoparquet/djibouti-latest_filter-road/raw/slice-0_edges.geoparquet",
            "results/geoparquet/djibouti-latest_filter-road/raw/slice-0_nodes.geoparquet",
        )
    )
