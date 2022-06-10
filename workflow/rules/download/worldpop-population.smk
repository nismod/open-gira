"""Download Worldpop population counts, constrained individual countries 2020 UN adjusted 
(100m resolution)

Reference
---------
https://www.worldpop.org/geodata/listing?id=79
"""

# the following worldpop.org request is sometimes responded to with 403 permission denied
# the default user-agent header of the requests library is 'python-requests/2.27.1'
# changing this to one from a browser will circumvent the access control and return 200
headers = {
    'user-agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:100.0) Gecko/20100101 Firefox/100.0',
}
r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m", headers=headers)
COUNTRY_CODES = [row["iso3"] for row in r.json()["data"]]

out_population = expand(
    os.path.join(config['output_dir'], "input", "population", "{country}_ppp_2020_UNadj_constrained.tif"),
    country=COUNTRY_CODES,
)


rule download_population:
    params:
        output_dir = config['output_dir'],
        code_country="{country}"
    output:
        os.path.join(config['output_dir'], "input", "population", "{country}_ppp_2020_UNadj_constrained.tif"),
    script:
        "../../scripts/download/scrape_url.py"
