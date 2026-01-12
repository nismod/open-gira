"""
Tropical cyclone tracks from Chia-Ying Lee at Columbia.

https://doi.org/10.1002/2017MS001186

Tracks available upon request as ragged netCDFs. These tracks have been
further processed with: https://github.com/thomas-fred/chaz-track-parser to
calibrate annual frequencies to historical record, per-basin. The calibrated
frequency tracks, in tabular format, are used as input here.
"""


rule parse_CHAZ:
    """
    Assume input is a symlink to most recent version of processed CHAZ on SoGE cluster.

    N.B. Only using Column Relative Humidity (CRH) (as opposed to Saturation Deficit) tracks.

    Test with:
    snakemake -c1 results/storm_tracks/CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050/0/tracks.geoparquet
    """
    input:
        # Zero-pad sample number to match input format
        # Parse scenario name, e.g.: CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050
        raw = lambda wildcards: (
            f"{wildcards.OUTPUT_DIR}/input/CHAZ/"
            f"CHAZ_genesis-CRH_SSP-{wildcards.CHAZ_SCENARIO.split('_SSP-')[1].split('_GCM')[0]}"
            f"_GCM-{wildcards.CHAZ_SCENARIO.split('_GCM-')[1].split('_epoch')[0]}"
            f"_epoch-{wildcards.CHAZ_SCENARIO.split('_epoch-')[1]}"
            f"_sample-{int(wildcards.SAMPLE):03d}.gpq"
        )
    output:
        processed = "{OUTPUT_DIR}/storm_tracks/{CHAZ_SCENARIO}/{SAMPLE}/tracks.geoparquet"
    script:
        "./parse_CHAZ.py"


rule slice_CHAZ:
    """
    Subset tracks to country convex hull plus a buffer.

    To test:
    snakemake -c1 results/power/by_country/PRI/storms/CHAZ-SSP-585-epoch-2050-GCM-UKESM1-0-LL/0/tracks.geoparquet
    """
    input:
        global_tracks="{OUTPUT_DIR}/storm_tracks/{CHAZ_SCENARIO}/{SAMPLE}/tracks.geoparquet",
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/{CHAZ_SCENARIO}/{SAMPLE}/tracks.geoparquet",
    resources:
        mem_mb=10000  # the global tracks file is fairly chunky
    script:
        "./slice_storm_tracks.py"

