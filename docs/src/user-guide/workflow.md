# Exploring output

The program will have created a great many files in subdirectories of `./results`.
Each of these will be the result of one of the steps in the workflow process:
1. Calculate a grid of bounding boxes for slicing the OpenStreetMap data
2. Filter the OpenStreetMap data to focus on major infrastructure components
3. Slice the OpenStreetMap data into smaller sections
4. Convert OpenStreetMap data to .geoparquet format
5. Add hazard information to infrastructure geoparquet
6. Join slices together to produce overall geoparquet file

These steps, along with the output produced at each stage, 
are described in the subsections of this chapter.
