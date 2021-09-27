configfile: 'config.yaml'

data_DIR = config['data_dir']
OUTPUT_DIR = config['output_dir']

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


rule filter_osm_data:
    input: os.path.join(DATA_DIR, '{pbf_file}.osm.pbf')
    output: os.path.join(DATA_DIR, '{pbf_file}-highway-core.osm.pbf')
    shell: 'osmium tags-filter {input} w/highway={filters} -o {output}'


rule convert_to_geoparquet:
    input:
        cmd='osm_to_pq.py',
        data=os.path.join(DATA_DIR, '{pbf_file}-highway-core.osm.pbf')
    output: os.path.join(DATA_DIR, '{pbf_file}-highway-core.geoparquet')
    shell: 'python {input.cmd} {input.data} {OUTPUT_DIR}'

rule clean:
    shell: 'rm -f data/*-highway-core.osm.pbf data/*.geoparquet'
