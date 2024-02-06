rule download_IBTrACS:
    output:
        "{OUTPUT_DIR}/input/IBTrACS/raw/v4.csv"
    shell:
        """
        wget --output-document {output} \
            https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.ALL.list.v04r00.csv
        """

"""
Test with:
snakemake -c1 results/input/IBTrACS/v4.csv
"""

rule parse_ibtracs:
    input:
        ibtracs_csv = "{OUTPUT_DIR}/input/IBTrACS/raw/v4.csv"
    output:
        ibtracs_parquet = "{OUTPUT_DIR}/storm_tracks/IBTrACS/0/tracks.geoparquet"
    script:
        "./parse_IBTrACS.py"

"""
To test:
snakemake -c1 results/storm_tracks/IBTrACS/0/tracks.geoparquet
"""


rule slice_ibtracs:
    input:
        global_tracks=rules.parse_ibtracs.output.ibtracs_parquet,
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/IBTrACS/0/tracks.geoparquet",
    script:
        "./slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/IBTrACS/0/tracks.geoparquet
"""
