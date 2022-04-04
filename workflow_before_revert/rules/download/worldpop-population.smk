"""Download Worldpop population counts, constrained individual countries 2020 UN adjusted 
(100m resolution)

Reference
---------
https://www.worldpop.org/geodata/listing?id=79
"""

r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m")
COUNTRY_CODES = [row["iso3"] for row in r.json()["data"]]

out_population = expand(
    os.path.join('data', "population", "{country}_ppp_2020_UNadj_constrained.tif"),
    country=COUNTRY_CODES,
)


rule download_population:
    input:
        out_population,


rule download_population_indiv:
    output:
        os.path.join('data', "population", "{code}_ppp_2020_UNadj_constrained.tif"),
    params:
        code_country = "{code}",
    script:
            os.path.join("..", "..", "scripts", "download", "scrape_url.py"
        )