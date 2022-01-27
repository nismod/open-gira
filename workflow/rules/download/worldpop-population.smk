"""Download Worldpop population counts, constrained individual countries 2020 UN adjusted 
(100m resolution)

Reference
---------
https://www.worldpop.org/geodata/listing?id=79
"""

out_population = expand(
    os.path.join(DATA_DIR, "population", "{country}_ppp_2020_UNadj_constrained.tif"),
    country=COUNTRY_CODES,
)


rule download_population:
    input:
        out_population

rule download_population_indiv:
    output:
        os.path.join(DATA_DIR, "population", "{code}_ppp_2020_UNadj_constrained.tif"),
    shell:
         "python3 "+os.path.join("workflow", "scripts", "download", "scrape_url.py")+" {wildcards.code}"


