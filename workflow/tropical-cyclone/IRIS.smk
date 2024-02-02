"""
Download the Imperial tropical cyclone storm tracks
"""

rule download_IRIS:
    """
    As of 20230626, this data is not publically available. You will need
    appropriate keys to access the files on the OUCE file store.
    """
    output:
        zip_file = "{OUTPUT_DIR}/input/IRIS/archive.zip"
    shell:
        """
        mkdir -p $(dirname {output.zip_file})
        scp /ouce-home/projects/mistral/iris/iris-data.zip {output.zip_file}
        """

"""
Test with:
snakemake -n -c1 -- results/input/IRIS/archive.zip
"""


rule extract_IRIS:
    input:
        zip_file = rules.download_IRIS.output.zip_file
    output:
        unzipped_dir = directory("{OUTPUT_DIR}/input/IRIS/iris-data/")
    shell:
        """
        unzip {input.zip_file} -d $(dirname {output.unzipped_dir})
        """

rule parse_IRIS:
    input:
        csv_dir="{OUTPUT_DIR}/input/IRIS/iris-data/event_sets/{IRIS_SCENARIO}/"
    output:
        parquet="{OUTPUT_DIR}/storm_tracks/IRIS-{IRIS_SCENARIO}/{SAMPLE}/tracks.geoparquet"
    script:
        "./parse_IRIS.py"

"""
Test with:
snakemake -c1 results/storm_tracks/IRIS_SSP1-2050/0/tracks.geoparquet
"""


rule slice_IRIS:
    input:
        global_tracks=rules.parse_IRIS.output.parquet,
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/IRIS-{IRIS_SCENARIO}/{SAMPLE}/tracks.geoparquet",
    resources:
        mem_mb=10000  # the global tracks file is fairly chunky
    script:
        "./slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/IRIS-PRESENT/0/tracks.geoparquet
"""
