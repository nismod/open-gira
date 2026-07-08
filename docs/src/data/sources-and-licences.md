# Data sources and licences

open-gira is built on open data, but "open" covers a range of licences, and the
obligations differ. This page is the register of the upstream datasets open-gira
can download, what each is used for, and the licence it comes under. It serves
two readers: someone deciding whether they may use an open-gira output for their
purpose, and the maintainers keeping the project's licensing defensible.

> **This page is an audit in progress.** The datasets below are enumerated from
> `config/config.yaml` and the workflow rules, and are therefore accurate as to
> *what open-gira uses*. The **licence** column, however, must be verified
> against each provider before it is relied upon — licences change, and a wrong
> entry here would be worse than an empty one. Entries marked *(verify)* have not
> yet been confirmed for this page.
>
> *(maintainer to-do: confirm each licence and its redistribution terms against
> the provider, add dataset version, approximate size, and typical download
> time, and resolve the flagged item below.)*

## Infrastructure and base data

| Dataset             | Used for                              | Source                                                        | Licence                    |
| ------------------- | ------------------------------------- | ------------------------------------------------------------- | -------------------------- |
| OpenStreetMap (via [Geofabrik](https://download.geofabrik.de/)) | Road and rail networks | Geofabrik extracts of OSM | ODbL 1.0                    |
| [gridfinder](https://github.com/carderne/gridfinder) | Electricity transmission/distribution | Predicted grid lines | *(verify)*                 |
| WRI Global Power Plant Database | Power plant locations         | World Resources Institute                                     | *(verify)*                 |
| [GADM](https://gadm.org/) 3.6 | Administrative boundaries      | gadm.org                                                       | **see flagged item below** |
| [Natural Earth](https://www.naturalearthdata.com/) | Context / coastlines | Natural Earth                                                 | Public domain              |
| GHSL (built-up, population) | Buildings exposure, population   | EC Joint Research Centre                                       | *(verify)*                 |

## Hazard data

| Dataset             | Hazard                | Source                                                        | Licence      |
| ------------------- | --------------------- | ------------------------------------------------------------- | ------------ |
| Aqueduct Flood      | River & coastal flood | WRI Aqueduct Flood Tool                                        | *(verify)*   |
| JRC river flood     | River flood           | EC Joint Research Centre                                       | *(verify)*   |
| Deltares coastal    | Coastal flood         | Deltares                                                       | *(verify)*   |
| IBTrACS             | Tropical cyclone (historic tracks) | NOAA                                             | *(verify)*   |
| STORM               | Tropical cyclone (synthetic)       | Synthetic track set                             | *(verify)*   |
| IRIS                | Tropical cyclone (synthetic)       | Synthetic track set                             | *(verify)*   |
| CHAZ                | Tropical cyclone (synthetic)       | Synthetic track set                             | *(verify)*   |
| Emanuel tracks      | Tropical cyclone (synthetic)       | Synthetic track set                             | *(verify)*   |
| Landslide (Arup)    | Landslide susceptibility           | Arup                                            | *(verify)*   |

The full, authoritative list of hazard file locations is in
`config/hazard_resource_locations/`; the source URLs in `config/config.yaml` are
the definitive record of what open-gira fetches.

## Flagged item: administrative boundaries

open-gira currently uses **GADM 3.6** for administrative boundaries. GADM's
licence permits academic and other non-commercial use but restricts
redistribution and commercial use. This has a practical consequence: open-gira
can *download* GADM for a user's own analysis, but **open-gira data releases that
embed GADM-derived boundaries may not be freely redistributable**.

Because a goal of open-gira is to publish reusable data products, this is a
constraint worth resolving rather than working around. One option under
consideration is [geoBoundaries](https://www.geoboundaries.org/), which is
released under terms friendlier to redistribution. Any change here should be made
deliberately, because it affects the identifiers and geometries downstream
analyses depend on.

*(maintainer to-do: decide on the boundary source for redistributable releases,
and record the decision and its rationale here.)*

## Using open-gira outputs

open-gira's own outputs are produced under the project's
[MIT licence](https://github.com/nismod/open-gira/blob/main/LICENSE), but an
output's *redistributability* is also governed by the licences of the upstream
data it derives from — most importantly OpenStreetMap's ODbL (which carries
share-alike and attribution obligations) and the boundary-data question above.
When you publish an open-gira output, check the obligations of every input that
went into it.

*(How-to guides for loading open-gira outputs in QGIS, GeoPandas, and CLIMADA
are planned; see the documentation outline.)*
