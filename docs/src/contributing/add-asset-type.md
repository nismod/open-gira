# Add an asset type and damage curve

Suppose open-gira builds a network containing a kind of asset it does not yet
estimate damages for — a new road class, a new rail structure, a substation. To
bring that asset into the damage calculation you need to do three things: declare
the asset type, give it a damage curve, and give it a rehabilitation cost. This
guide walks through each, using the existing road and rail assets as the model.

The pieces fit together like this: a network feature is tagged with an
`asset_type`; the damage step looks up a **damage curve** (hazard intensity →
damage fraction) for that asset type; and the cost step multiplies the damage
fraction by a **rehabilitation cost** to get a monetary loss. If any of the three
is missing for an asset type you request, the workflow will tell you.

## 1. Declare the asset type

Asset types are defined in `src/open_gira/assets.py` as string-valued enums.
Each sector has a subclass of `Assets`; for example, road assets:

```python
class RoadAssets(Assets):
    """Typology of road assets"""

    BRIDGE = "road_bridge"
    MOTORWAY = "road_motorway"
    TRUNK = "road_trunk"
    PRIMARY = "road_primary"
    # ...
    UNPAVED = "road_unpaved"
```

To add an asset type, add a member to the appropriate subclass (or create a new
`Assets` subclass for a new sector). The **value** is the `asset_type` string
that will appear in networks and in file names, so follow the existing
convention of `<sector>_<name>`:

```python
    CYCLEWAY = "road_cycleway"
```

`Assets.implemented_assets()` collects the values from every subclass, and the
configuration check in `workflow/Snakefile` uses it to validate the
`direct_damages.asset_types` you request — so once declared, the new type is
immediately recognised.

Your network-creation logic must actually tag features with the new
`asset_type` for it to be used; where that tagging happens depends on the sector
(see the network-creation code for that sector).

## 2. Add a damage curve

A damage curve is a CSV under `config/damage_curves/<hazard_type>/<asset_type>.csv`
— for example `config/damage_curves/flood/road_cycleway.csv`. The file format,
documented in `config/damage_curves/flood/README.md`, is:

1. **Source header**, one or more lines each prefixed with `#`, citing where the
   curve comes from. This is not optional politeness — it is how the provenance
   of a damage estimate stays attached to it.
2. **Column headings** on the first non-comment line. For flooding this is
   exactly:
   ```
   inundation_depth_(m),damage_fraction
   ```
3. **The data**, one record per line: hazard intensity, then damage fraction in
   the range 0–1.

For example:

```
# Source: Author et al., 2024, Table 2
# https://doi.org/...
inundation_depth_(m),damage_fraction
0.0,0.0
0.5,0.2
1.0,0.33
2.0,0.5
```

The intensity units and the heading depend on the hazard type; match the other
curves in the same `config/damage_curves/<hazard_type>/` directory.

## 3. Add a rehabilitation cost

Rehabilitation costs live in `config/rehab_costs/<sector>.csv`, one row per
asset type, giving the cost to rebuild per unit (for linear assets, per
kilometre — the field name used internally is `rehab_cost_USD_per_km`). Add a row
for the new asset type, again with a `#`-prefixed source header describing where
the figure comes from, following the format of the existing
`config/rehab_costs/road.csv` and `rail.csv`.

## 4. Request it

With the asset type declared, its curve in place, and its cost recorded, add it
to `direct_damages.asset_types` in `config/config.yaml` (or leave that list
empty to use every implemented asset type), and request a damages target as in
the [transport/flooding tutorial](../tutorials/wales-roads-flooding.md). The new
asset type will now contribute to the damage and cost estimates.

## 5. Add a test

Following the [testing conventions](index.md#how-the-tests-are-structured), add
or extend a test so the new asset type's damage calculation is exercised on the
sample dataset. This keeps the new curve from silently breaking later.
