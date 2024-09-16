# Multi-modal network integration

Combine different transport network modes (road, rail, maritime) together into a
connected global network, centered on a single country.

## Description

1. We require input land networks to connect. These can be created with other workflows in open-gira or you may bring your own.

Using open-gira you can create a road network for Thailand with the following:
```bash
snakemake --cores all results/thailand-latest_filter-road-primary/edges.gpq
```
See more details in [road.md](./road.md) and [rail.md](./rail.md).

Or make a composite network from
```bash
snakemake --cores all results/composite_network/south-east-asia-road/edges.gpq
```
For more details, see [composite-networks.md](./composite-networks.md).

However they're created, input land networks should be placed in the following locations:
```python
road_network_nodes = "{OUTPUT_DIR}/input/networks/road/{PROJECT_SLUG}/nodes.gpq",
road_network_edges = "{OUTPUT_DIR}/input/networks/road/{PROJECT_SLUG}/edges.gpq",
rail_network_nodes = "{OUTPUT_DIR}/input/networks/rail/{PROJECT_SLUG}/nodes.gpq",
rail_network_edges = "{OUTPUT_DIR}/input/networks/rail/{PROJECT_SLUG}/edges.gpq",
```

1. We also require input maritime network data:

Network nodes, including ports for import, export and trans-shipment.
```python
nodes = "{OUTPUT_DIR}/input/networks/maritime/nodes.gpq",
```
Sample node data:
```
     id infra                                   name iso3 Continent_Code                   geometry
  port3  port                Aberdeen_United Kingdom  GBR             EU  POINT (-2.07662 57.13993)
 port84  port               Avonmouth_United Kingdom  GBR             EU  POINT (-2.70877 51.49812)
port121  port                   Barry_United Kingdom  GBR             EU  POINT (-3.24712 51.40467)
```

Network edges with costs (used for routing over). The referenced port nodes are
appended with `in` or `out` as the network is bi-directional.
```python
edges_no_geom = "{OUTPUT_DIR}/input/networks/maritime/edges_by_cargo/maritime_base_network_general_cargo.pq",
```
Sample edge data:
```
     from_id         to_id from_infra to_infra      mode from_iso3 to_iso3   distance_km  cost_USD_t_km from_port   to_port
   port0_out   port1099_in       port     port  maritime       AUS     ZAF  14296.116346       0.002385     port0  port1099
port1001_out   port1099_in       port     port  maritime       BRA     ZAF   9016.208416       0.003095  port1001  port1099
port1005_out   port1099_in       port     port  maritime       ITA     ZAF  10973.998701       0.004198  port1005  port1099
```

Network edges without costs, but with a geometry that can be used for visualising flows over.
```python
edges_visualisation = "{OUTPUT_DIR}/input/networks/maritime/edges.gpq",
```
Sample visualisation edge data:
```
    from_id         to_id from_infra  to_infra               id                                           geometry
  maritime0  maritime1265   maritime  maritime  maritimeroute_0  LINESTRING (11.00000 -19.00002, 5.51647 -19.58...
maritime213  maritime1265   maritime  maritime  maritimeroute_1  LINESTRING (9.99999 -30.00001, 4.79550 -25.083...
maritime964  maritime1265   maritime  maritime  maritimeroute_2  LINESTRING (-10.00000 -10.00002, -7.58272 -12....
```

## Configuration

- Review and amend `config/config.yaml`:
  - `study_country_iso_a3`: the ISO 3 letter code of the country in question
  - `road_cost_USD_t_km`: cost in USD per tonne per km of road freight
  - `road_cost_USD_t_h`: cost in USD per tonne per hour of road freight
  - `road_default_speed_limit_km_h`: speed limit in km per hour if no other is available
  - `rail_cost_USD_t_km`: cost in USD per tonne per km of rail freight
  - `rail_cost_USD_t_h`: cost in USD per tonne per hour of rail freight
  - `rail_average_freight_speed_km_h`: speed assumption for rail freight in km per hour
  - `intermodal_cost_USD_t`: cost in USD per tonne of changing from one network mode to another (e.g. road to rail)

## Creation

Here's an example command to create a multi-modal network, bringing together land networks under `project-thailand` and joining them to a maritime network:
```bash
snakemake --cores all results/multi-modal_network/project-thailand/edges.gpq
```