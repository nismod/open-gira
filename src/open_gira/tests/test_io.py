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


class TestPackingFloatsAsInts:
    """
    Use xarray to round trip floats via netCDF int16.
    """

    def test_round_trip_to_netcdf(self):

        # maximum error of 0.3%
        relative_tolerance = 0.003

        arrays_to_test = [
            np.array([-23.32, -2.2, 1, 24.21, 128.3, 243.1245]),
            np.linspace(0.1, 81.234, 1500),
            np.linspace(0.138, 978.234, 500),
            np.linspace(-999.34, 10012.21, 40),
            np.random.rand(100),
            np.exp(np.random.rand(100)),
            np.logspace(-2, 20, 10, base=2) + np.random.rand(10),
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

                try:
                    assert np.allclose(data, from_disk.data.values, rtol=relative_tolerance, atol=0, equal_nan=True)
                except AssertionError:
                    ratio = from_disk.data.values / data
                    error = np.abs(ratio - 1)
                    raise AssertionError(f"{np.nanmax(error)=:.5f} exceeds {relative_tolerance=:.5f}")
