# STORM IBTrACS data

The synthetic cyclone storm data is retrieved from https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085 and is used to simulate storms on the power network. The corresponding return periods are obtained from https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164 which are currently only used to generate the grid resolution of the intersection analysis. More details on the generated data can be found in the paper "Generation of a global synthetic tropical cyclone hazard dataset using STORM" by N. Bloemendaal, I. D. Haigh, H. de Moel, S. Muis, R. J. Haarsma & C. J. H. Jeroen from 2020 (https://www.nature.com/articles/s41597-020-0381-2).

The data is split into 6 basins:
 - EP: East Pacific
 - NA: North Atlantic
 - NI: North Indian
 - SI: South Indian
 - SP: South Pacific
 - WP: West Pacific

The two letter abbreviation(s) is/are selected in the `config.yaml` file for analysis. The aforementioned paper also visually presents the storm basis.