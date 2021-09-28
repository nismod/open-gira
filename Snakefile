configfile: 'config.yaml'

DATA_DIR = config['data_dir']
OUTPUT_DIR = config['output_dir']
AQUEDUCT_DIR = config['aqueduct_dir']

links = [
    "motorway",
    "motorway_link",
    "trunk",
    "trunk_link",
    "primary_link",
    "secondary",
    "secondary_link"
]
filters = ','.join(links)

FULL_PBF_FILE = os.path.join(DATA_DIR, "{slug}.osm.pbf")
PBF_FILE = os.path.join(DATA_DIR, "{slug}-highway-core.osm.pbf")
GEOPARQUET_FILE = PBF_FILE.replace(".osm.pbf", ".geoparquet")

GEOPARQUET_SPLITS_FILE = GEOPARQUET_FILE.replace(
    ".geoparquet", "_splits.geoparquet"
).replace(DATA_DIR, OUTPUT_DIR)

PARQUET_SPLITS_FILE = GEOPARQUET_SPLITS_FILE.replace(".geoparquet", ".parquet")


rule filter_osm_data:
    shell: 'osmium tags-filter {input} w/highway={filters} -o {output}'
    input:
        FULL_PBF_FILE,
    output:
        PBF_FILE,


rule convert_to_geoparquet:
    input:
        cmd='osm_to_pq.py',
    shell: 'python {input.cmd} {input.data} {DATA_DIR}'
        data=PBF_FILE,
    output:
        GEOPARQUET_FILE,


rule network_hazard_intersection:
    input:
        cmd='network_hazard_intersection.py',
        csv=os.path.join(AQUEDUCT_DIR, 'aqueduct_river.csv')
        network=GEOPARQUET_FILE,
    output:
    shell: 'python {input.cmd} {input.network} {AQUEDUCT_DIR} {OUTPUT_DIR}'
        geoparquet=GEOPARQUET_SPLITS_FILE,
        parquet=PARQUET_SPLITS_FILE,


rule clean:
    shell: 'rm -rf data/*-highway-core.osm.pbf data/*.geoparquet outputs/'
