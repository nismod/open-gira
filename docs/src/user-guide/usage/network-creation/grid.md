# Electricity grid

To create an electricity grid we rely heavily on [gridfinder](https://gridfinder.rdrn.me/)
data. This dataset provides transmission and distribution edges with
substantial global coverage. It also contains a set of 'targets' or electricity
consuming areas, derived from [Night Time Lights](https://www.earthdata.nasa.gov/learn/backgrounders/nighttime-lights)
(NTL) satellite imagery. Our other major data source for electricity grid creation is the World Resources Institute's (WRI)
[powerplants database](https://datasets.wri.org/dataset/globalpowerplantdatabase).

The workflow currently divides network creation by country. One may request one
country, or several. Note that neighbouring countries' networks are _not_
connected with one another.

## Description

- Download gridfinder electricity grid data (line edges, 'target' consuming polygons)
- Create global 'targets' file, where each electricity consuming target is annotated with the population and GDP of the enclosed area.
- Download WRI powerplant data
- Download GDP per capita data
- Download population data
- Construct an electricity grid for each country to be studied where nodes are power generators or target consumers. Power consumption is proportional to the GDP enclosed by each target polygon.

## Configuration

There aren't currently any options to configure when creating an electricity
grid.

## Running

Here's an example grid creation command for Puerto Rico:
```bash
snakemake --cores all -- results/power/by_country/PRI/network/edges.geoparquet
```

Note that the folder name under `by_county` is an [ISO 3166 Alpha-3 country
code](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3#Officially_assigned_code_elements),
specifying the country in question.
