import os
from urllib.parse import urlparse

import requests
import sys

# Get the page with the table
r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m")


for ii,row in enumerate(r.json()['data']):
    country_ident = row['iso3']
    country_url = "https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m?iso3="+country_ident
    country_data = requests.get(country_url).json()
    link = country_data['data'][0]['files'][0]

    try:
        output_path = snakemake.output  # get snakemake output naming
        if output_path[ii].lower() != (os.path.join("data","population",os.path.basename(urlparse(link).path))).lower():
            raise RuntimeError("Naming convention for population files is different")  # check out_population naming
    except NameError:
        # If "snakemake" doesn't exist then must be running from the command line.
        raise RuntimeError("Run from snakemake file, else you will need to specify approx 200 files in cmd line")

    fname = output_path[ii]  # successful naming (checks passed)

    if os.path.exists(fname):  # check if file already exists
        continue
    print(fname)
    with open(fname, 'wb') as fh:
        data_request = requests.get(link, stream=True)
        for chunk in data_request.iter_content(chunk_size=1024*1024*4):
            fh.write(chunk)
