"""Download Worldpop population counts, constrained individual countries 2020 UN adjusted 
(100m resolution)

Reference
---------
https://www.worldpop.org/geodata/listing?id=79
"""

r = requests.get(
    "https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m",
    headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"},
)
COUNTRY_CODES = [row["iso3"] for row in r.json()["data"]]

out_population = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "population",
        "{country}_ppp_2020_UNadj_constrained.tif",
    ),
    country=COUNTRY_CODES,
)


rule download_population:
    params:
        output_dir=config["output_dir"],
        code_country="{country}",
    output:
        os.path.join(
            config["output_dir"],
            "input",
            "population",
            "{country}_ppp_2020_UNadj_constrained.tif",
        ),
    script:
        "../../scripts/download/scrape_url.py"


rule download_population_all:
    input:
        out_population,
