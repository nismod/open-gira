"""
Tropical cyclone tracks from Kerry Emanuel at MIT / WindRiskTech LLC.

Tracks available upon request as ragged (padded) Matlab. These tracks have been
further processed with: https://github.com/thomas-fred/emanuel-track-parser to
calibrate annual frequencies to historical record, per-basin. The calibrated
frequency tracks, in tabular format, are used as input here.
"""


rule parse_emanuel:
    """
    Read and parse track set.

    N.B. ~200 years total, so less than one full millennium sample.

    Test with:
    snakemake -c1 results/storm_tracks/emanuel_ssp-585_gcm-ukmo6_epoch-2050/0/tracks.geoparquet
    """
    input:
        raw = "{OUTPUT_DIR}/input/emanuel/{EMANUEL_SCENARIO}_tracks.gpq"
    output:
        processed = "{OUTPUT_DIR}/storm_tracks/{EMANUEL_SCENARIO}/0/tracks.geoparquet"
    script:
        "./parse_emanuel.py"


rule slice_emanuel:
    """
    Subset tracks to country convex hull plus a buffer.

    To test:
    snakemake -c1 results/power/by_country/PRI/storms/emanuel_ssp-585_gcm-ukmo6_epoch-2050/0/tracks.geoparquet
    """
    input:
        global_tracks="{OUTPUT_DIR}/storm_tracks/{EMANUEL_SCENARIO}/{SAMPLE}/tracks.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{EMANUEL_SCENARIO}/{SAMPLE}/tracks.geoparquet",
    resources:
        mem_mb=10000
    script:
        "./slice_storm_tracks.py"

