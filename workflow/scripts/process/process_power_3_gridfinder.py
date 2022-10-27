"""Indexes all gridfinder values"""
import geopandas


if __name__ == "__main__":
    output_dir = snakemake.config["output_dir"]  # type: ignore
    gridfinder_path = snakemake.input.gridfinder  # type: ignore
    global_boxes_path = snakemake.input.global_boxes  # type: ignore
    output_path = snakemake.output.gridfinder  # type: ignore
    box_id = snakemake.wildcards.BOX  # type: ignore

    boxes = geopandas.read_file(global_boxes_path) \
        .set_index("box_id")
    box = boxes.loc[[f"box_{box_id}"], :]
    gridfinder = geopandas.read_file(gridfinder_path).reset_index(names="source_id")
    gridfinder_box = gridfinder.sjoin(box).rename(columns={"index_right": "box_id"})
    gridfinder_box.to_parquet(output_path)
