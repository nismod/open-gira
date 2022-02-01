import os
from urllib.parse import urlparse

import requests
import sys


country_ident = sys.argv[1]
country_url = (
    "https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m?iso3=" + country_ident
)
country_data = requests.get(country_url).json()
link = country_data["data"][0]["files"][0]

fname = os.path.join(
    "data", "population", f"{country_ident}_ppp_2020_UNadj_constrained.tif"
)

if not os.path.exists(fname):  # check if file already exists
    print(fname)
    with open(fname, "wb") as fh:
        data_request = requests.get(link, stream=True)
        for chunk in data_request.iter_content(chunk_size=1024 * 1024 * 4):
            fh.write(chunk)
