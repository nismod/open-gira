"""Download Worldpop population counts, constrained individual countries 2020 UN adjusted 
(100m resolution)

Reference
---------
https://www.worldpop.org/geodata/listing?id=79
"""

out_population = expand(
    os.path.join(DATA_DIR, "population", "{country}_ppp_2020_UNadj_constrained.tif"), 
    country=COUNTRY_CODES)

rule download_population:
    output: out_population
    script: "../scripts/scrape_url.py"
