configfile: 'config.yaml'

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
    input: 'data/tanzania-latest.osm.pbf',
    output: 'data/tanzania-latest-highway-core.osm.pbf'
    shell: 'osmium tags-filter {input} w/highway={filters} -o {output}'
        
