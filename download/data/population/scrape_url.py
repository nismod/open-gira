import os
from urllib.parse import urlparse

import requests
from bs4 import BeautifulSoup

# Get the page with the table
# TODO parameterise DOI or dataset id
r = requests.get("https://www.worldpop.org/ajax/geolisting/doi?doi=10.5258/SOTON/WP00685")



for row in r.json():
    data_page = requests.get("https://www.worldpop.org/geodata/summary", params={'id':row['id']})
    soup = BeautifulSoup(data_page.text, 'html.parser')
    links = soup.find(id='files')
    link = links.find('a').get('href')
    fname = os.path.basename(urlparse(link).path)
    if os.path.exists(fname):
        continue
    print(fname)
    with open(fname, 'wb') as fh:
        data_request = requests.get(link, stream=True)
        for chunk in data_request.iter_content(chunk_size=1024*1024*4):
            fh.write(chunk)
