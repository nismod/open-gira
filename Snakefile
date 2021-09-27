configfile: 'config.yaml'

INPUT_DIR = config['input_dir']
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
    input: os.path.join(INPUT_DIR, '{pbf_file}.osm.pbf')
    output: os.path.join(OUTPUT_DIR, '{pbf_file}-highway-core.osm.pbf')    
    shell: 'osmium tags-filter {input} w/highway={filters} -o {output}'

rule clean:
    shell: 'rm data/*-highway-core.osm.pbf'
        
