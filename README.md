# OpenStreetMap / Aqueduct Flood exposure

Use `nismod/snail` to demonstrate infrastructure network intersections with hazards as an
exposure calculation.

Goals: 
- automated pipeline for reproducible exposure analysis anywhere in the world.
- map of a large area showing exposure at different return periods
- charts/stats of exposure per admin region (by road/rail) per hazard type, scenario, epoch

## Data

### Gridfinder

```
pip install zenodo_get
mkdir -p data/gridfinder
cd data/gridfinder
zenodo_get -d 10.5281/zenodo.3628142
```

### WRI Aqueduct

Use `nismod/aqueduct`. 

TODO update to provide snail raster dataset metadata.

### OpenStreetMap

Testing with downloads from http://download.geofabrik.de/asia.html - Laos and Bangladesh.

TODO consider planet and diffs.

