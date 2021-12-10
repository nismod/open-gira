# Download data files
CYCLONE_REGIONS = ["EP","NA","NI","SI","SP","WP"]
r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m")
COUNTRY_CODES = [row['iso3'] for row in r.json()['data']]

out_fixed = expand(os.path.join(DATA_DIR, "stormtracks", "fixed", "STORM_FIXED_{param}_{let}.nc"), let = CYCLONE_REGIONS, param = ["RETURN_PERIODS","TC_WIND_SPEEDS"])
out_events = expand(os.path.join(DATA_DIR, "stormtracks", "events", "STORM_DATA_IBTRACS_{let}_1000_YEARS_{num}.txt"), let = CYCLONE_REGIONS, num=list(range(0,10)))
out_population = expand(os.path.join(DATA_DIR, "population", "{country}_ppp_2020_UNadj_constrained.tif"), country = COUNTRY_CODES)
out_GDP = expand(
    os.path.join(DATA_DIR, "GDP", "{filename}"), 
    filename=[
        "admin_areas_GDP_HDI.nc", 
        "GDP_per_capita_PPP_1990_2015_v2.nc",
        "GDP_PPP_30arcsec_v3.nc",
        "GDP_PPP_1990_2015_5arcmin_v2.nc",
        "HDI_1990_2015_v2.nc",
        "pedigree_GDP_per_capita_PPP_1990_2015_v2.nc",
        "pedigree_HDI_1990_2015_v2.nc"
    ])
out_powerplant = os.path.join(DATA_DIR, "powerplants", "global_power_plant_database.csv")
out_gridfinder = os.path.join(DATA_DIR, "gridfinder", "grid.gpkg")
out_adminboundaries = os.path.join(DATA_DIR, "adminboundaries", "gadm36.gpkg")
out_adminboundaries_codes = expand(
    os.path.join(DATA_DIR, "adminboundaries", "gadm36_{code}.gpkg"),
    code = COUNTRY_CODES)

rule download_all:
    input: 
        out_fixed,
        out_events,
        out_population,
        out_GDP,
        out_powerplant,
        out_gridfinder,
        out_adminboundaries_codes  # choose between whole world or per-country

rule download_stormtracks_fixed:
    output: out_fixed
    shell: 
        """
        wget \
            -P data/stormtracks/fixed \
            -N \
            -i workflow/scripts/storm_fixed_return.txt \
            --no-check-certificate
        """

rule download_stormtracks_events:
    output: out_events
    shell: 
        """
        wget \
            -P data/stormtracks/events \
            -N \
            -i workflow/scripts/storm_tracks.txt \
            --no-check-certificate \
            --content-disposition
        unzip -o data/stormtracks/events/STORM_DATA3.zip -d data/stormtracks/events
        """

rule download_population:
    output: out_population
    script: "../scripts/scrape_url.py"

rule download_GDP:
    output: out_GDP
    shell: 
        """
        mkdir -p data/GDP \
        cd data/GDP \
        wget https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.dk1j0/download --output-document=doi_10.5061_dryad.dk1j0__v2.zip &&
        unzip -o doi_10.5061_dryad.dk1j0__v2.zip
        """

rule download_powerplants:
    output: out_powerplant
    shell: 
        """
        mkdir -p data/powerplants
        cd data/powerplants
        wget https://wri-dataportal-prod.s3.amazonaws.com/manual/global_power_plant_database_v_1_3.zip --output-document=global_power_plant_database_v_1_3.zip
        unzip -o global_power_plant_database_v_1_3.zip
        """

# download gridfinder data
rule download_gridfinder:
    output: out_gridfinder
    shell: 
        """
        mkdir -p data/gridfinder
        cd data/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """

rule download_gadm:
    output: out_adminboundaries
    shell: 
        """
        mkdir -p data/adminboundaries
        cd data/adminboundaries
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip \
            --output-document=gadm36_gpkg.zip
        unzip -o gadm36_gpkg.zip
        """

# download admin boundaries (per country)
rule download_gadm_by_country:
    output: os.path.join(DATA_DIR, "adminboundaries", "gadm36_{code}.gpkg"),
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_{code}_gpkg.zip \
            --output-document=data/adminboundaries/gadm36_{code}_gpkg.zip
        unzip -o data/adminboundaries/gadm36_{code}_gpkg.zip -d data/adminboundaries
        """
