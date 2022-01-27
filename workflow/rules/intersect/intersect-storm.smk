"""Extract hazard data and intersect. Examine statistics

"""
import os


region_box = expand(os.path.join('data', 'intersection', 'regions', '{region}_boxes.txt'), region=REGIONS)
region_grid = expand(os.path.join('data', 'intersection', 'regions', '{region}_grid.gpkg'), region=REGIONS)
TC_years = expand(os.path.join("data","intersection", "storm_data", "all_winds", "log", "__winds_completed_r{region}_s{sample}_y{year}.txt"), region=REGIONS, sample=SAMPLES, year=YEARS)
stat_csv = os.path.join("data", "intersection", "combined_storm_statistics.csv")


rule intersect_regions_indiv:
    input:
        os.path.join('data', 'processed', 'world_boxes_metadata.txt'),
        #os.path.join("data", "processed", "world_boxes.gpkg"),  # removed because opening on QGIS tampers with metadata
        os.path.join("data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"),
        expand(os.path.join("data", "processed", "all_boxes", "{box_id}", "gridfinder_{box_id}.gpkg"), box_id=all_boxes)
    output:
        os.path.join('data', 'intersection', 'regions', '{region}_boxes.txt'),
    shell:
        "python3 "+os.path.join("workflow", "scripts", "intersect", "intersect_1_regions.py")+" {wildcards.region}"


rule intersect_regions:
    input:
        region_box


rule intersect_grid_indiv:
    input:
        expand(os.path.join("data", "processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"), box_id=all_boxes),
        os.path.join("data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"),
        os.path.join('data', 'intersection', 'regions', '{region}_boxes.txt')
    output:
        os.path.join("data", "intersection", "regions", "{region}_grid.gpkg")
    shell:
        "python3 "+os.path.join("workflow", "scripts", "intersect", "intersect_2_gridmaker.py")+" {wildcards.region}"


rule intersect_grid:
    input:
        region_grid


rule intersect_winds_indiv:
    input:
        os.path.join('data', 'processed', 'world_boxes_metadata.txt'),
        os.path.join("data", "intersection", "regions", "{region}_grid.gpkg"),
        os.path.join("data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"),
        os.path.join("data", "stormtracks", "events", "STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt")
    output:
        [os.path.join("data","intersection", "storm_data", "all_winds", "log", "__winds_completed_r{region}_s{sample}"+f"_y{year}.txt") for year in YEARS]  # use as pseudo file for csv (otherwise snakemake deletes output files everytime)
    shell:
        "python3 "+os.path.join("workflow", "scripts", "intersect", "intersect_3_winds.py")+" {wildcards.region} {wildcards.sample} "+'"""'+str(YEARS)+'"""'


rule intersect_wind:
    input:
        TC_years


rule intersection_gdploss:
    input:
        os.path.join("data","intersection", "storm_data", "all_winds", "log", "__winds_completed_r{region}_s{sample}_y{year}.txt"),
        os.path.join("data", "intersection", "regions", "{region}_grid.gpkg"),
        [os.path.join("data","processed", "all_boxes", f"{box_id}", f"network_with_gdp_{box_id}.gpkg") for box_id in all_boxes],
        [os.path.join("data","processed", "all_boxes", f"{box_id}", f"edge_gdp_sorted_{box_id}.txt") for box_id in all_boxes],
        [os.path.join("data","processed", "all_boxes", f"{box_id}", f"targets_{box_id}.gpkg") for box_id in all_boxes]
    output:
        os.path.join("data", "intersection","storm_data", "damages", "storm_r{region}_s{sample}_y{year}.txt")
    shell:
        "python3 "+os.path.join("workflow", "scripts", "intersect", "intersect_4_gdploss.py")+" {wildcards.region} {wildcards.sample} {wildcards.year} {operationfind}"


rule merge_all_stats:
    """Use this rule for the combined stats file"""
    input:
         expand(os.path.join("data", "intersection","storm_data", "damages", "storm_r{region}_s{sample}_y{year}.txt"), region=REGIONS, sample=SAMPLES, year=YEARS),
    output:
        stat_csv,
    shell:
        "python3 " +  os.path.join("workflow", "scripts", "intersect", "stat_merger.py")
