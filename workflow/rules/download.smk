# download storm tracks data
rule download_stormtracks_fixed:
	output: out_fixed
	shell: "wget -P data/stormtracks/fixed -N -i data/stormtracks/storm_fixed_return.txt --no-check-certificate"


rule download_stormtracks_events:
	output: out_events
	shell: """wget -P data/stormtracks/events -N -i data/stormtracks/storm_tracks.txt --no-check-certificate --content-disposition &&
			  unzip -o data/stormtracks/events/STORM_DATA3.zip -d data/stormtracks/events"""
			  

# download population data
rule download_population:
	output: out_population
	script: "../scripts/scrape_url.py"
	

# download GDP data
rule download_GDP:
	output: out_GDP
	shell: """mkdir -p data/GDP &&
			  cd data/GDP &&
			  wget https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.dk1j0/download --output-document=doi_10.5061_dryad.dk1j0__v2.zip &&
			  unzip -o doi_10.5061_dryad.dk1j0__v2.zip"""
			  
	
# download powerplants data
rule download_powerplants:
	output: out_powerplant
	shell: """mkdir -p data/powerplants &&
			  cd data/powerplants &&
			  wget https://wri-dataportal-prod.s3.amazonaws.com/manual/global_power_plant_database_v_1_3.zip --output-document=global_power_plant_database_v_1_3.zip &&
			  unzip -o global_power_plant_database_v_1_3.zip"""
	

# download gridfinder data
rule download_gridfinder:
	output: out_gridfinder
	shell: """mkdir -p data/gridfinder &&
			  cd data/gridfinder &&
			  zenodo_get 10.5281/zenodo.3628142"""
			  
			 
# download admin boundaries
rule download_adminboundaries:
	output: out_adminboundaries
	shell: """mkdir -p data/adminboundaries &&
			  cd data/adminboundaries &&
			  wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip --output-document=gadm36_gpkg.zip &&
			  unzip -o gadm36_gpkg.zip"""	
			  

# download admin boundaries (per country)
rule download_adminboundaries_codes:
	output: out_adminboundaries_codes
	run: 
		for code in COUNTRY_CODES:
			shell("""wget https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_{code}_gpkg.zip --output-document=data/adminboundaries/gadm36_{code}_gpkg.zip &&
			  unzip -o data/adminboundaries/gadm36_{code}_gpkg.zip -d data/adminboundaries""")
	
	
			  
