import os
import tempfile

import numpy as np
import xarray as xr

from open_gira.io import netcdf_packing_parameters


class TestScaleFactorAndOffset:
    """
    Test the packing of floats into integers on disk.
    """

    def test_regression(self):
        scale_factor, offset, fill_value = netcdf_packing_parameters(0, 80, 16)

        np.allclose(
            np.array([scale_factor, offset, fill_value]),
            np.array([0.0012207, 40.000610, -2**15])
        )


class TestPackingFloatsAsInts:
    """
    Use xarray to round trip floats via netCDF int16.
    """

    def test_round_trip_to_netcdf(self):

        # maximum error of 0.3%
        relative_tolerance = 0.003

        arrays_to_test = [
            np.array([-23.32, -2.2, 1, np.nan, 24.21, 69.25293649, 128.3, 243.1245]),
            np.array([np.nan, np.nan, 1, np.nan, 24.21, np.nan, 128.3, np.nan]),
            np.linspace(0.138, 78.234, 500),
            np.linspace(0.1, 81.234, 1500),
            np.linspace(-999.34, 10012.21, 40),
            np.random.rand(100),
            np.exp(np.random.rand(100)),
            np.logspace(-2, 7, 10, base=2) + 1/7
        ]

        for data in arrays_to_test:

            da = xr.DataArray(
                data,
                name="data",
                coords={"x": range(len(data))},
                dims=["x"]
            )
            scale_factor, offset, fill_value = netcdf_packing_parameters(np.nanmin(data), np.nanmax(data), 16)

            with tempfile.TemporaryDirectory() as temp_dir:
                serialised_path = os.path.join(temp_dir, "serialised.nc")

                da.to_netcdf(
                    serialised_path,
                    encoding={
                        "data": {
                            "dtype": "int16",
                            "scale_factor": scale_factor,
                            "add_offset": offset,
                            "_FillValue": fill_value
                        }
                    }
                )

                from_disk = xr.load_dataset(serialised_path)
                try:
                    assert np.allclose(data, from_disk.data.values, rtol=relative_tolerance, atol=0, equal_nan=True)
                except AssertionError:
                    ratio = from_disk.data.values / data
                    relative_error = np.abs(ratio - 1)
                    max_relative_error = np.nanmax(relative_error)
                    raise AssertionError(f"{max_relative_error=:.5f} exceeds {relative_tolerance=:.5f}")
