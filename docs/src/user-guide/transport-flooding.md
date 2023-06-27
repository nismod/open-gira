# Transport - flooding

The steps in the workflow process:
1. Download OpenStreetMap data
2. Filter OpenStreetMap data to focus on major infrastructure components
3. Determine the bounding box from the OpenStreetMap data
4. Download hazard raster files
5. Clip hazard raster files to bounding boxes determined in step 3
6. Calculate a grid of bounding boxes for slicing the OpenStreetMap data
7. Slice the OpenStreetMap data into smaller sections
8. Convert OpenStreetMap data to .geoparquet format
9. Add hazard information to infrastructure geoparquet
10. Join slices together to produce overall geoparquet file

These steps, along with the output produced at each stage,
are described in the subsections of this chapter.

These steps are summarised in the digital acyclic graph for `slice_count: 1`, for just the
`tanzania-latest` infrastructure and `aqueduct-coast` hazard data.

[![DAG of the workflow for
the Tanzania dataset and coast flooding data](./img/DAG-simple.png)](./img/DAG-simple.png)
