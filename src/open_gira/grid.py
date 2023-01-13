import pandas as pd


def weighted_allocation(
    nodes: pd.DataFrame,
    *,
    variable_col: str,
    weight_col: str,
    component_col: str,
    asset_col: str,
    source_name: str,
    sink_name: str,
) -> pd.DataFrame:
    """
    For each network component, allocate total variable capacity of sources
    to sinks (consuming nodes), weighted by some feature of the sinks. The sum
    of the variable for all sources and sinks in a component should equal zero.

    Args:
        nodes: Should contain, `variable_col`, `weight_col`, `asset_col` and `component_col`
        variable_col: Name of column in `nodes` to distribute from `source` to `sink`
        weight_col: Name of column in `nodes` to weight allocation to sinks by
        component_col: Name of column in `nodes` to group sources and sinks by
        asset_col: Name of column in `nodes` to differentiate sources and sinks by
        source_name: Categorical for sources in `nodes[asset_col]`
        sink_name: Categorical for sinks in `nodes[asset_col]`

    Returns:
        Sink nodes with variable allocated from source nodes by weight
    """

    # find the sum of variable for each component
    c_variable_sum = nodes.loc[
        nodes[asset_col] == source_name,
        [variable_col, component_col]
    ].groupby(component_col).sum().reset_index()

    # subset to sinks
    sinks = nodes[nodes[asset_col] == sink_name]

    # find the sum of weights for each component
    c_weight_sum = sinks.loc[:, [weight_col, component_col]].groupby(component_col).sum().reset_index()

    # merge in the component sums for variable and weight
    c_variable_col = f"_component_{variable_col}"
    sinks = sinks.merge(
        c_variable_sum.rename(columns={variable_col: c_variable_col}),
        how="left",
        on=component_col
    )
    c_weight_col = f"_component_{weight_col}"
    sinks = sinks.merge(
        c_weight_sum.rename(columns={weight_col: c_weight_col}),
        how="left",
        on=component_col
    )

    # ensure every sink has a numeric (non-NaN) entry for the component sum of variable
    sinks[c_variable_col] = sinks[c_variable_col].fillna(0)

    # reallocate variable to sinks, by weight within components
    sinks.loc[:, variable_col] = -1 * sinks[c_variable_col] * sinks.loc[:, weight_col] / sinks[c_weight_col]

    return sinks
