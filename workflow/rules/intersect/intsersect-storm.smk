"""Extract hazard data and intersect. Examine statistics

"""
import os
import pandas as pd

def examine_hurr_list(region, sample):
    """Returns a list of all hurricane unique identifiers #_# where the first # is the sample and the second # is the number of the hurricane in that sample"""
    filename = os.path.join(DATA_DIR, "stormtracks", "events", f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt")
    df = pd.read_csv(filename, header=None)
    df.columns = ['year', 'month', 'number', 'step','basin','lat','lon','pressure','wind','radius','cat','landfall','dis_land']
    nums = list(df['number'].unique())
    return [f"{sample}_{int(num)}" for num in nums]


RP_nums = list(range(10,100,10))+list(range(100, 1000, 100))+list(range(1000,11000,1000))  # required for nc to tif conversion
RP_stats = ['conf_5', 'conf_95', 'mean', 'stdev']
REGIONS = ["WP"]  # must ensure this covers the COUNTRY_CODES
SAMPLES = [0]  # Range from 0 to 10
EXAMINE_HURR = ["0_0", "0_1", "0_2", "0_3", "0_4"]  # Hurricanes to examine  # issues for EXAMINE_HURR if more than 1 region (since different number of hurricanes very possible)
EXAMINE_HURR = examine_hurr_list("WP", 0)[:30]  # TODO 30 for testing
print(EXAMINE_HURR)

opfind = True  # find operational fraction of targets, requires extra approx 40% time

RP_tifs = expand(os.path.join(DATA_DIR, "stormtracks", "fixed", "extracted",  "STORM_FIXED_RETURN_PERIODS_{region}_rp{RP_num}_{RP_stat}.tif"), RP_num = RP_nums, RP_stat = RP_stats, region = REGIONS)
WS_csvs = expand(os.path.join("data","intersection", 'TC_c{code}_r{region}_s{sample}.csv'), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES)
stat_csv = os.path.join("data", "intersection", "combined_storm_statistics.csv")

rule intersect_all:
    input:
        RP_tifs,
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}", "world_edges_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}","world_region_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data","intersection", "storm_data", "storm_{nh}","world_targets__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
        expand(os.path.join("data", "intersection", "storm_data", "storm_{nh}", "storm_c{code}_r{region}_s{sample}_n{nh}.txt"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),


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
        #os.path.join("data","intersection", 'storms_all', '{nh}__c{code}_r{region}_s{sample}.txt'),
    shell:
        "python3 " +  os.path.join(WORKFLOW_DIR, "scripts", "intersect", "wind_speeds.py") + " {wildcards.code} {wildcards.region} {wildcards.sample}"


rule hurr_examine:
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


rule merge_all_stats:
    """Use this rule if the combined stats file needs update or is incomplete. The other rules should ensure this does not happen though."""
    input:
        expand(os.path.join("data", "intersection", "storm_data", "storm_{nh}", "storm_c{code}_r{region}_s{sample}_n{nh}.txt"), code=COUNTRY_CODES, region=REGIONS, sample=SAMPLES, nh=EXAMINE_HURR),
    output:
        stat_csv,
    shell:
        "python3 " +  os.path.join(WORKFLOW_DIR, "scripts", "intersect", "stat_merger.py")
