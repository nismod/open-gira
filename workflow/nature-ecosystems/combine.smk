"""Extract NbS opportunities as polygon patches, intersected with Hydrobasins
and attributed with costs and benefits.

Then join with infrastructure/hazard EAD for relevant combinations:
- mangrove / coastal flooding
- slope vegetation / landslide
- river basin afforestation / river flooding


`data.extracts[0].bbox` is xmin, ymin, xmax, ymax in degrees:
```json
{"directory": ".", "extracts": [{"bbox": [-180.0,-60.0,-172.0,-57.06666666666667],"output": "slice-0.osm.pbf"}]}
```
"""
import geopandas

def fname_sector(fname: str) -> str:
    if "filter-rail" in fname:
        sector = "rail"
    elif "filter-road" in fname:
        sector = "road"
    else:
        assert False, f"Unexpected sector in {fname}"
    return sector

def read_ead(fname: str, sector: str) -> geopandas.GeoDataFrame:
    ead = geopandas.read_parquet(fname)
    ead_cols = [c for c in ead.columns if "_EAD" in c]
    asset_cols = ["geometry"]
    rename_cols = {col: f"{col}__{sector}" for col in ead_cols}
    ead = ead[asset_cols + ead_cols].rename(columns=rename_cols)
    return ead

rule slice_raster:
    input:
        tiff="{OUTPUT_DIR}/input/nbs-suitability/{DATASET}/{KEY}.tif",
        json="{OUTPUT_DIR}/json/{DATASET}_extracts/{SLICE_SLUG}.geojson",
    output:
        tiff="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/{KEY}.tif",
    shell:
        """
        set -ex

        mkdir -p $(dirname {output.tiff})

        # pull out bounding box coords into bash array
        # bbox is [xmin, ymin, xmax, ymax]
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json}  | sed 's/"//g'))"

        # then use gdal_translate with projwin to avoid warping
        # projwin is [ulx, uly, lrx, lry]
        gdal_translate \
            -projwin ${{COORDS[0]}} ${{COORDS[3]}} ${{COORDS[2]}} ${{COORDS[1]}} \
            -co compress=lzw \
            --config GDAL_CACHEMAX 50% \
            {input.tiff} {output.tiff}

        # gdalwarp gives identical grids, but hard to guarantee pixel alignment with source data
        # -te is [xmin, ymin, xmax, ymax]
        # gdalwarp \
        #     -te ${{COORDS[@]}} \
        #     -co COMPRESS=LZW \
        #     --config GDAL_CACHEMAX 50% \
        #     {input.tiff} {output.tiff}
        """

rule slice_stack:
    input:
        biodiversity_benefit_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_BioBenefit_9s.tif",
        carbon_benefit_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_CarbonBenefit_9s.tif",
        planting_cost_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_PlantingCost_9s.tif",
        regeneration_cost_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_RegenCost_9s.tif",
        mangrove_planting_cost_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/ManPlantCost_9s.tif",
        mangrove_regeneration_cost_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/ManRegenCost_9s.tif",
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

rule mangrove_with_ead:
    input:
        opportunities="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/mangrove.parquet",
        risks=expand(
            "{{OUTPUT_DIR}}/direct_damages/{{DATASET}}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD_and_cost_per_RP/{{SLICE_SLUG}}.geoparquet",
            HAZARD_SLUG=[
                "hazard-aqueduct-coast",
                "hazard-coastal-deltares",
            ],
            FILTER_SLUG=[
                "filter-rail",
                "filter-road-tertiary",
            ],
        ),
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/mangrove_with_EAD.parquet",
    run:
        import geopandas
        from open_gira.io import write_empty_frames

        MANGROVE_EFFECT_M = 5000
        options = geopandas.read_parquet(input.opportunities)

        if options.empty:
            write_empty_frames(output.parquet)
            return


        buf_geom = options.geometry.to_crs("ESRI:54009").buffer(MANGROVE_EFFECT_M).to_crs("EPSG:4326")
        buf = geopandas.GeoDataFrame(
            data = options.reset_index()[["feature_id"]],
            geometry =  buf_geom.reset_index(drop=True)
        )

        for fname in input.risks:
            sector = fname_sector(fname)
            ead = read_ead(fname, sector)
            joined = buf.sjoin(ead, how="left").drop(columns=["geometry", "index_right"])
            options_damages = joined.groupby("feature_id").sum()
            options = options.join(options_damages)
        options.to_parquet(output.parquet)

rule landslide_opportunities:
    input:
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        nbs_stack="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.parquet",
        landslide_slope_potential_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_LandslideNbS_123_9s.tif",
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/landslide_slope_vegetation.parquet",
    script:
        "./combine.landslide_opportunities.py"

rule landslide_with_ead:
    input:
        opportunities="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/landslide_slope_vegetation.parquet",
        risks=expand(
            "{{OUTPUT_DIR}}/direct_damages/{{DATASET}}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD_and_cost_per_trigger/{{SLICE_SLUG}}.geoparquet",
            HAZARD_SLUG=[
                "hazard-landslide-arup",
            ],
            FILTER_SLUG=[
                "filter-rail",
                "filter-road-tertiary",
            ],
        ),
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/landslide_slope_vegetation_with_EAD.parquet",
    run:
        import geopandas
        from open_gira.io import write_empty_frames
        LANDSLIDE_EFFECT_M = 3000
        options = geopandas.read_parquet(input.opportunities)

        if options.empty:
            write_empty_frames(output.parquet)
            return

        buf_geom = options.geometry.to_crs("ESRI:54009").buffer(LANDSLIDE_EFFECT_M).to_crs("EPSG:4326")
        buf = geopandas.GeoDataFrame(
            data = options.reset_index()[["feature_id"]],
            geometry =  buf_geom.reset_index(drop=True)
        )

        for fname in input.risks:
            sector = fname_sector(fname)
            ead = read_ead(fname, sector)
            joined = buf.sjoin(ead, how="left").drop(columns=["geometry", "index_right"])
            options_damages = joined.groupby("feature_id").sum()
            options = options.join(options_damages)
        options.to_parquet(output.parquet)

rule river_opportunities:
    input:
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        nbs_stack="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/nbs_stack.parquet",
        tree_potential_tif="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/G_PotentialNonCoastalTreeNBS_9s.tif",
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/river_basin_afforestation.parquet",
    script:
        "./combine.river_opportunities.py"

rule river_with_ead:
    input:
        opportunities="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/river_basin_afforestation.parquet",
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        risks=expand(
            "{{OUTPUT_DIR}}/direct_damages/{{DATASET}}_{FILTER_SLUG}/{HAZARD_SLUG}/EAD_and_cost_per_RP/{{SLICE_SLUG}}.geoparquet",
            HAZARD_SLUG=[
                "hazard-aqueduct-river",
                "hazard-river-jrc",
            ],
            FILTER_SLUG=[
                "filter-rail",
                "filter-road-tertiary",
            ],
        ),
    output:
        parquet="{OUTPUT_DIR}/slices/{DATASET}_nbs/{SLICE_SLUG}/river_basin_afforestation_with_EAD.parquet",
    run:
        import geopandas
        from open_gira.io import write_empty_frames
        options = geopandas.read_parquet(input.opportunities)

        if options.empty:
            write_empty_frames(output.parquet)
            return

        # union options so there's one multipolygon per basin - mean per_ha values
        options = options.dissolve(
            by=['HYBAS_ID', 'GID_0', 'GID_1', 'GID_2'],
            aggfunc={
                'area_m2': 'sum',
                'area_ha': 'sum',
                'biodiversity_benefit': 'mean',
                'carbon_benefit_t_per_ha': 'mean',
                'planting_cost_usd_per_ha': 'mean',
                'regen_cost_usd_per_ha': 'mean',
            },
            dropna=False,
        ).reset_index(level=['GID_0','GID_1','GID_2'])

        # in future, could read with bbox option, if written with write_covering_bbox=True
        # for now, read all and select with bounds using spatial index
        basins = geopandas.read_parquet(input.hydrobasins_adm, columns=["HYBAS_ID", "geometry"]).reset_index()
        xmin, ymin, xmax, ymax = options.total_bounds
        basins = basins.cx[xmin:xmax, ymin:ymax]

        # associate all damages within a basin to full option area within the basin
        for fname in input.risks:
            sector = fname_sector(fname)
            ead = read_ead(fname, sector)
            joined = basins.sjoin(ead, how="left").drop(columns=["geometry", "index_right"])
            options_damages = joined.groupby("HYBAS_ID").sum()
            options = options.join(options_damages)
        options.to_parquet(output.parquet)


rule combine_with_ead:
    input:
        parquet=expand(
            "{{OUTPUT_DIR}}/slices/{{DATASET}}_nbs/slice-{SLICE_SLUG}/{{NBS}}_with_EAD.parquet",
            SLICE_SLUG=range(config["slice_count"])
        )
    output:
        parquet="{OUTPUT_DIR}/nbs-adaptation/{DATASET}/{NBS}_with_EAD.parquet",
    run:
        # read from possibly-mixed schemas - would be nice to handle these as a single parquet directory
        # with common schema, but would need to set this up at write time, and pass to `write_empty_files`
        pds = pq.ParquetDataset(input.parquet)
        nonempty = []
        for pf_name in pds.files:
            pf = pq.ParquetFile(pf_name)
            if "HYBAS_ID" in pf.schema.names:
                nonempty.append(pf_name)
        pt = pq.read_table(nonempty)

        # Write to combined geoparquet file
        df = pt.to_pandas()
        df.geometry = geopandas.GeoSeries.from_wkb(df.geometry)
        gdf = geopandas.GeoDataFrame(df).set_crs(epsg=4326)
        gdf.to_parquet(output.parquet)

rule combine_all_with_ead:
    input:
        parquet=expand(
            "{{OUTPUT_DIR}}/nbs-adaptation/{{DATASET}}/{NBS}_with_EAD.parquet",
            NBS=["river_basin_afforestation", "mangrove", "landslide_slope_vegetation"]
        )
    output:
        touch("{OUTPUT_DIR}/nbs-adaptation/{DATASET}/combined.flag")
