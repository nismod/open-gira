# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files
rule network_hazard_intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}.geoparquet",
        tifs=glob(os.path.join(config['hazard_data_dir'], "*.tif")),
        hazard_csv=config["hazard_csv"],
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.parquet",
    script:
        "../scripts/network_hazard_intersection.py"

rule test_network_hazard_intersection:
    input:
        geoparquet=expand(
            os.path.join(
                config['output_dir'],
                'splits',
                f"{config['dataset']}_filter-{filter_slug}_slice-{{i}}_hazard-{hazard_slug}.geoparquet"
            ),
            i=range(config['slice_count'])
        ),
        parquet=expand(
            os.path.join(
                config['output_dir'],
                'splits',
                f"{config['dataset']}_filter-{filter_slug}_slice-{{i}}_hazard-{hazard_slug}.parquet"
            ),
            i=range(config['slice_count'])
        )