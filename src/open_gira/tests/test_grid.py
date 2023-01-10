import pandas as pd
import numpy as np

from open_gira.grid import allocate_power_to_targets


class Test_allocate_power_to_targets:

    def test_allocate_power_to_targets(self):
        df = pd.DataFrame(
            index=range(5),
            data={
                "component_id": [1, 1, 1, 2, 2],
                "asset_type": ["source", "target", "target", "target", "source"],
                "power_mw": [20, np.nan, np.nan, np.nan, 10],
                "gdp": [np.nan, 3, 1, 1, np.nan],
            }
        )

        generated = allocate_power_to_targets(df, "gdp")

        assert all(pd.Series([20, -15, -5, -10, 10]) == generated.power_mw)
