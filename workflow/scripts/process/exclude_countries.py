"""Writes a list of countries which are in the gadm36.gpkg file but do not have population .tifs to file
"""

import json
import os
import glob
import fiona

try:
    output_dir = snakemake.params['output_dir']
except:
    output_dir = sys.argv[1]



with fiona.open(
    os.path.join(output_dir, "input", "adminboundaries", "gadm36_levels.gpkg"), "r", layer=0
) as src_code:
    gadm36_countries = []
    for feature in src_code:
        gadm36_countries.append(feature["properties"]["GID_0"])


print("finding tif codes")
files = glob.glob(os.path.join(output_dir, "input", "population", "*.tif"))
pop_countries = [file[file.find("_ppp") - 3 : file.find("_ppp")] for file in files]

exclude_country_list = []

print("Finding excluded countries")
for country in gadm36_countries:
    if country not in pop_countries:
        exclude_country_list.append(country)

print("writing to file")
folder_path = os.path.join(output_dir, "power_processed")
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

with open(os.path.join(folder_path, "exclude_countries.txt"), "w") as file:
    json.dump(exclude_country_list, file)
print("Exported excluded countries")
