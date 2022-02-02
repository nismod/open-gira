# Exploring output

The program will have created a great many files in subdirectories of `./results`.
Each of these will be the result of one of the steps in the workflow process:
1. Slice the OpenStreetMap data into smaller sections
2. Filter the OpenStreetMap data to focus on major infrastructure components
3. Convert OpenStreetMap data to .geoparquet format
4. Add hazard information to infrastructure geoparquet
5. Join slices together to produce overall geoparquet file

These steps, along with the output produced at each stage, 
are described in the subsections of this chapter.
