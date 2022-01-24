"""Extract hazard data and intersect. Examine statistics

"""

# 
# rule region_boxes:
#     input:
#         os.path.join('data', 'processed', 'world_boxes_metadata.txt'),
#         os.path.join("data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc")
#     output:
#         os.path.join('data', 'processed', 'regions', '{region}_boxes.txt')
#     shell:
#         "python3 intersect_1_regions.py "+str(wildcards.region)
# 
# 
# rule wind_speeds_extract:
#     input:
#         os.path.join('data', 'processed', 'regions', '{region}_boxes.txt')
#     output:
#         os.path.join("data","intersection", "storm_data", "all_winds", 'TC_r{region}_s{sample}_n{nh_}.csv')