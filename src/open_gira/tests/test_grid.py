import pandas as pd
import numpy as np

from open_gira.grid import weighted_allocation


class Test_allocate_power_to_targets:
    def test_weighted_allocation(self):
        df = pd.DataFrame(
            index=range(5),
            data={
                "component_id": [1, 1, 1, 2, 2],
                "asset_type": ["source", "target", "target", "target", "source"],
                "power_mw": [20, np.nan, np.nan, np.nan, 10],
                "gdp": [np.nan, 3, 1, 1, np.nan],
            },
        )

        generated = weighted_allocation(
            df,
            variable_col="power_mw",
            weight_col="gdp",
            component_col="component_id",
            asset_col="asset_type",
            source_name="source",
            sink_name="target",
        )

        assert all(pd.Series([-15, -5, -10]) == generated["power_mw"])
        assert all(pd.Series([20, 20, 10]) == generated["_component_power_mw"])
        assert all(pd.Series([4, 4, 1]) == generated["_component_gdp"])
