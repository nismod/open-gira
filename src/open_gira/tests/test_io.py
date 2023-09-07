import os
import tempfile

import numpy as np
import xarray as xr

from open_gira.io import scale_factor_and_offset


class TestScaleFactorAndOffset:
    """
    Test the packing of floats into integers on disk.
    """

    def test_regression(self):
        scale_factor, offset = scale_factor_and_offset(0, 80, 16)

        np.allclose(
            np.array([scale_factor, offset]),
            np.array([0.0012207, 40.000610])
        )

    def test_round_trip_to_netcdf(self):
        """
        Use xarray to round trip via packed int16.
        """

        arrays_to_test = [
            np.array([-23.32, -2.2, 0, 24.21, 128.3, 27000]),
            np.array([-12.2365324, -2.1212, 0, 42.123, 128.6543, 27000.2363]),
            np.random.rand(100),
            np.random.rand(100) * 100_000,
            np.logspace(-2, 4, 10),
            np.logspace(-2, 9, 10, base=2),
            np.linspace(-1000, 10000, 40),
        ]

        for data in arrays_to_test:

            da = xr.DataArray(data, name="data", coords={"x": range(len(data))})
            scale_factor, offset = scale_factor_and_offset(data.min(), data.max(), 16)

            with tempfile.TemporaryDirectory() as temp_dir:
                serialised_path = os.path.join(temp_dir, "serialised.nc")

                da.to_netcdf(
                    serialised_path,
                    encoding={"data": {"dtype": "int16", "scale_factor": scale_factor, "add_offset": offset}}
                )

                from_disk = xr.load_dataset(serialised_path)

                np.allclose(
                    data,
                    from_disk.data.values
                )
