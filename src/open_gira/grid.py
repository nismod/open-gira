import geopandas as gpd


def allocate_power_to_targets(nodes: gpd.GeoDataFrame, weighting: str) -> gpd.GeoDataFrame:
    """
    For each network component, allocate total generating capacity of sources
    to targets (consuming nodes), weighted by some feature of the targets'.

    Here's some returned `nodes` weighted by `gdp`. This function has mutated
    the `power_mw` cells for target `asset_types`:
                        id   source_id asset_type  component_id          gdp  power_mw
    37      source_37_1030  WRI1028006     source             6          NaN  6.000000
    956    target_914_1030         NaN     target             6 8.256908e+07 -4.848396
    993    target_951_1030         NaN     target             6 1.441001e+06 -0.084615
    994    target_952_1030         NaN     target             6 2.087217e+06 -0.122560
    1000   target_958_1030         NaN     target             6 2.353413e+06 -0.138191
    1001   target_959_1030         NaN     target             6 1.342598e+07 -0.788364
    1009   target_967_1030         NaN     target             6 3.044258e+05 -0.017876

    Args:
        nodes: Should contain `id`, `component_id`, `power_mw`, `asset_type`
            and column referenced by `weighting`
        weighting: Name of column in `nodes` to weight targets by

    Returns:
        Nodes with generating capacity allocated to consuming nodes by weighting
    """

    power = "power_mw"

    # for each component, allocate generation weighted by GDP of targets
    for c_id in set(nodes.component_id):

        c_mask: pd.Series = nodes.component_id == c_id

        c_total_capacity: float = nodes.loc[(nodes.asset_type == "source") & c_mask, power].sum()

        c_target_mask: pd.Series = (nodes.asset_type == "target") & c_mask

        c_total_weight: float = nodes.loc[c_target_mask, weighting].sum()

        nodes.loc[c_target_mask, power] = \
            -1 * c_total_capacity * (nodes.loc[c_target_mask, weighting] / c_total_weight)

    return nodes
