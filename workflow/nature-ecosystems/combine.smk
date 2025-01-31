"""Extract NbS opportunities as polygon patches, intersected with Hydrobasins
and attributed with costs and benefits.

TODO join with infrastructure/hazard EAD for relevant combinations:
- mangrove / coastal flooding
- slope vegetation / landslide
- river basin afforestation / river flooding


`data.extracts[0].bbox` is xmin, ymin, xmax, ymax in degrees:
```json
{"directory": ".", "extracts": [{"bbox": [-180.0,-60.0,-172.0,-57.06666666666667],"output": "slice-0.osm.pbf"}]}
```
"""


rule slice_raster:
    input:
        tiff="{OUTPUT_DIR}/input/nbs-suitability/{DATASET}/{KEY}.tif",
        json="{OUTPUT_DIR}/json/{DATASET}_extracts/{SLICE_SLUG}.geojson",
    output:
        tiff="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/{KEY}.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.tiff})

        # pull out bounding box coords into bash array
        # bbox is [xmin, ymin, xmax, ymax]
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json}  | sed 's/"//g'))"

        # then use gdal_translate with projwin to avoid warping
        # projwin is [ulx, uly, lrx, lry]
        gdal_translate \
            -projwin ${{COORDS[0]}} ${{COORDS[3]}} ${{COORDS[2]}} ${{COORDS[1]}} \
            -co compress=lzw \
            --config GDAL_CACHEMAX=50% \
            {input.tiff} {output.tiff}

        # gdalwarp gives identical grids, but hard to guarantee pixel alignment with source data
        # -te is [xmin, ymin, xmax, ymax]
        # gdalwarp \
        #     -te ${{COORDS[@]}} \
        #     -co COMPRESS=LZW \
        #     --config GDAL_CACHEMAX=50% \
        #     {input.tiff} {output.tiff}
        """

rule slice_stack:
    input:
        biodiversity_benefit_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_BioBenefit_9s.tif",
        carbon_benefit_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_CarbonBenefit_9s.tif",
        planting_cost_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_PlantingCost_9s.tif",
        regeneration_cost_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_RegenCost_9s.tif",
    output:
        zarr=directory("{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.zarr"),
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.parquet",
    script:
        "./combine.slice_stack.py"

rule mangrove_opportunities:
    input:
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        nbs_stack="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.parquet",
        mangrove_potential_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/ManRestorClass_9s.tif",
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/mangrove.parquet",
    script:
        "./combine.mangrove_opportunities.py"

rule landslide_opportunities:
    input:
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        nbs_stack="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.parquet",
        landslide_slope_potential_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_LandslideNbS_123_9s.tif",
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/landslide_slope_vegetation.parquet",
    script:
        "./combine.landslide_opportunities.py"

rule river_opportunities:
    input:
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        nbs_stack="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.parquet",
        tree_potential_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_PotentialNonCoastalTreeNBS_9s.tif",
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/river_basin_afforestation.parquet",
    script:
        "./combine.river_opportunities.py"
