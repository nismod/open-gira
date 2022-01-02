"""Extract hazard data and intersect. Examine statistics

"""
import os
import pandas as pd
import glob as gb

def examine_hurr_list(region, sample):
    """Returns a list of all hurricane unique identifiers #_# where the first # is the sample and the second # is the number of the hurricane in that sample"""
    filename = os.path.join(DATA_DIR, "stormtracks", "events", f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt")
    df = pd.read_csv(filename, header=None)
    df.columns = ['year', 'month', 'number', 'step','basin','lat','lon','pressure','wind','radius','cat','landfall','dis_land']
    nums = list(df['number'].unique())
    return [f"{sample}_{int(num)}" for num in nums]

def nh_req(nh_lst_all):
    """This function ensures that uncompleted intersections are done only. Important to have this function otherwise parallel computions wont work."""
    nh_lst_rem = nh_lst_all.copy()
    print("Removing pre-run storms")
    for stormnh in nh_lst_all:
        storm_folder = os.path.join("data", "intersection","storm_data", f"storm_{stormnh}")
        storm_files = gb.glob(os.path.join(storm_folder, "*.txt"))

        if len(storm_files) == 1:
            nh_lst_rem.remove(stormnh)
            print(f"{stormnh} already intersected")
    return nh_lst_rem  # all remaining


RP_nums = list(range(10,100,10))+list(range(100, 1000, 100))+list(range(1000,11000,1000))  # required for nc to tif conversion
RP_stats = ['conf_5', 'conf_95', 'mean', 'stdev']
REGIONS = ["WP"]  # must ensure this covers the COUNTRY_CODES
SAMPLES = [0]  # Range from 0 to 10
EXAMINE_HURR = ["0_0", "0_1", "0_2", "0_3", "0_4"]  # Hurricanes to examine  # issues for EXAMINE_HURR if more than 1 region (since different number of hurricanes very possible)
EXAMINE_HURR = examine_hurr_list("WP", 0)[:10]  # TODO # for testing
EXAMINE_HURR_req = nh_req(EXAMINE_HURR)
print(EXAMINE_HURR_req)

opfind = True  # find operational fraction of targets, requires extra approx 40% time

RP_tifs = expand(os.path.join(DATA_DIR, "stormtracks", "fixed", "extracted",  "STORM_FIXED_RETURN_PERIODS_{region}_rp{RP_num}_{RP_stat}.tif"), RP_num = RP_nums, RP_stat = RP_stats, region = REGIONS)
WS_csvs = expand(os.path.join("data","intersection", 'TC_c{code}_r{region}_s{sample}.csv'), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES)
stat_csv = os.path.join("data", "intersection", "combined_storm_statistics.csv")


rule intersect_all:
    """Use this when not many storm files need to be intersected. It has to load up the processed data for each storm intersection"""
    input:
        RP_tifs,
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_edges_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}","world_region_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}","world_targets__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data", "intersection", "storm_data", "storm_{nh}", "storm_c{code}_r{region}_s{sample}_n{nh}.txt"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),


rule intersect_all_parallel:
    """uses hurr_examine_parallel to load the intersecting data once to then paraellise (within the .py file) the intersection but still loads the data for every storm. Possible issues when using multiple cores already and then .py wanting to use more, see: cpu_count()-2"""
    input:
        RP_tifs,
        expand(os.path.join("data", "intersection","storm_data", "storm_{nh}", "storm_c{code}_r{region}_s{sample}_n{nh}_p.txt"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_edges_affected__storm_c{code}_r{region}_s{sample}_n{nh}_p.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_region_affected__storm_c{code}_r{region}_s{sample}_n{nh}_p.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_targets__storm_c{code}_r{region}_s{sample}_n{nh}_p.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        stat_csv,

rule nc_to_tif:  # TODO update so that output here is an input in another .py (will be hazard_intersect_sortedgdp.py once return periods added)
    input:
        os.path.join(DATA_DIR, "stormtracks", "fixed",  "STORM_FIXED_RETURN_PERIODS_{region}.nc"),
    output:
        os.path.join(DATA_DIR, "stormtracks", "fixed", "extracted",  "STORM_FIXED_RETURN_PERIODS_{region}_rp{RP_num}_{RP_stat}.tif"),
    shell:
        "python3 " +  os.path.join(DATA_DIR, "workflow", "scripts", "intersect", "stormtracks_nc_to_tif.py") + " {wildcards.region}"


rule wind_speeds_extract:
    input:
        os.path.join("data","adminboundaries", "gadm36_{code}.gpkg"),
    output:
        #os.path.join("data","intersection", "grid_{code}.gpkg"),
        os.path.join("data","intersection", 'TC_c{code}_r{region}_s{sample}.csv'),
    shell:
        "python3 " +  os.path.join(WORKFLOW_DIR, "scripts", "intersect", "wind_speeds.py") + " {wildcards.code} {wildcards.region} {wildcards.sample}"


rule hurr_examine:
    """Examines each storm individually, ideal if loading files are small"""
    input:
        #os.path.join("data","intersection", "grid_{code}.gpkg"),
        WS_csvs,
        os.path.join("data","processed", "world_network_with_gdp.gpkg"),
        os.path.join("data","processed", "world_targets.gpkg"),
        os.path.join("data","intersection", "TC_c{code}_r{region}_s{sample}.csv"),
        os.path.join("data","processed","edge_gdp_sorted.txt"),
    output:
        os.path.join("data", "intersection","storm_data", "storm_{nh}", "storm_c{code}_r{region}_s{sample}_n{nh}.txt"),
        os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_edges_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"),
        os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_region_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"),
        os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_targets__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"),
    shell:
        "python3 " +  os.path.join(WORKFLOW_DIR, "scripts", "intersect", "hazard_intersect_sortedgdp.py") + " {wildcards.code} {wildcards.region} {wildcards.sample} {wildcards.nh} " + str(opfind)


rule hurr_examine_parallel:
    """Examines all storm inputs (in parallel from .py), ideal if loading files are large"""
    input:
        #os.path.join("data","intersection", "grid_{code}.gpkg"),
        WS_csvs,
        os.path.join("data","processed", "world_network_with_gdp.gpkg"),
        os.path.join("data","processed", "world_targets.gpkg"),
        os.path.join("data","processed","edge_gdp_sorted.txt"),
    output:
        [os.path.join("data", "intersection","storm_data", f"storm_{s_idx}", "storm_c{code}_r{region}_s{sample}_n"+f"{s_idx}_p.txt") for s_idx in EXAMINE_HURR_req],
        [os.path.join("data","intersection", "storm_data", f"storm_{s_idx}", "world_edges_affected__storm_c{code}_r{region}_s{sample}_n"+f"{s_idx}_p.gpkg") for s_idx in EXAMINE_HURR_req],
        [os.path.join("data","intersection", "storm_data", f"storm_{s_idx}", "world_region_affected__storm_c{code}_r{region}_s{sample}_n"+f"{s_idx}_p.gpkg") for s_idx in EXAMINE_HURR_req],
        [os.path.join("data","intersection", "storm_data", f"storm_{s_idx}", "world_targets__storm_c{code}_r{region}_s{sample}_n"+f"{s_idx}_p.gpkg") for s_idx in EXAMINE_HURR_req],
    shell:
        "python3 " +  os.path.join(WORKFLOW_DIR, "scripts", "intersect", "hazard_intersect_sortedgdp_parallel2.py") + " {wildcards.code} {wildcards.region} {wildcards.sample} "+' """' +str(EXAMINE_HURR_req)+'""" '  + str(opfind)


rule merge_all_stats:
    """Use this rule for the combined stats file (for parallel _p files)"""
    input:
        expand(os.path.join("data", "intersection", "storm_data", "storm_{nh}", "storm_c{code}_r{region}_s{sample}_n{nh}_p.txt"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
    output:
        stat_csv,
    shell:
        "python3 " +  os.path.join(WORKFLOW_DIR, "scripts", "intersect", "stat_merger.py")
