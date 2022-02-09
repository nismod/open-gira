# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files
def txt_target_files(wildcards):
    """
    Return a list of the basenames of the download targets in the file_list.txt files for each hazard directory
    """
    with open(os.path.join(f"{config['output_dir']}", "input", wildcards.HAZARD_SLUG, "file_list.txt"), "r") as txt:
        out = [os.path.basename(x) for x in txt.readlines()]
    return [x.replace('\n', '') for x in out]

rule intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}.geoparquet",
        tifs=lambda wildcards: expand(
            os.path.join(
                f"{config['output_dir']}",
                "input",
                f"{wildcards.HAZARD_SLUG}",
                dataset_slug,
                f"{{filename}}"
            ),
            filename=txt_target_files(wildcards=wildcards)
        )
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.parquet",
    script:
        "../scripts/intersection.py"

# Test via target file:
# snakemake --cores all results/splits/tanzania-latest_filter-highway-core_slice-0_hazard-aqueduct-coast.geoparquet
