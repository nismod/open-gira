import os
import logging
import tempfile

import numpy as np
import pytest
import xarray as xr

from open_gira.io import netcdf_packing_parameters, \
    bit_pack_dataarray_encoding, bit_pack_dataset_encoding


LOGGER = logging.getLogger(__name__)


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

    def test_all_zeros(self):
        """
        Data is all zeros, so min and max are the same. Should return no transformation (1, 0).
        """
        scale_factor, offset, fill_value = netcdf_packing_parameters(0, 0, 16)

        np.allclose(
            np.array([scale_factor, offset, fill_value]),
            np.array([1, 0, -2**15])
        )

    def test_all_fill_value(self):
        """
        No variance in the data, and all the values are the typical fill_value. How awkward!
        Should return a fill_value at the other end of the integer range.
        """
        scale_factor, offset, fill_value = netcdf_packing_parameters(-2**15, -2**15, 16)

        np.allclose(
            np.array([scale_factor, offset, fill_value]),
            np.array([1, 0, 2**15])
        )


class TestPackingFloatsAsInts:
    """
    Use xarray to round trip floats via netCDF int16.
    """

    # maximum error of 0.3%
    relative_tolerance = 0.003

    np.random.seed(0)

    @staticmethod
    def round_trip_within_tolerance(da: xr.DataArray, encoding: dict, relative_tolerance: float) -> bool:
        """
        Given a DataArray, write it to disk with a stretch and scale integer
        encoding scheme, read it back and compare to the copy in-memory.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            serialised_path = os.path.join(temp_dir, "serialised.nc")
            da.to_netcdf(serialised_path, encoding=encoding)

            from_disk = xr.load_dataset(serialised_path)
            result = np.allclose(da.data, from_disk.data.values, rtol=relative_tolerance, atol=0, equal_nan=True)

        if result:
            return True
        else:
            ratio = from_disk.data.values / da.data
            relative_error = np.abs(ratio - 1)
            max_relative_error = np.nanmax(relative_error)
            print(f"{max_relative_error=:.5f} exceeds {relative_tolerance=:.5f}")
            return False

    def test_netcdf_packing_parameters_round_trip(self):
        """
        Round trip values using the core function, netcdf_packing_parameters directly.
        """
        arrays_to_test = [
            np.array([0]),
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
            encoding={
                "data": {
                    "dtype": "int16",
                    "scale_factor": scale_factor,
                    "add_offset": offset,
                    "_FillValue": fill_value
                }
            }
            assert self.round_trip_within_tolerance(da, encoding, self.relative_tolerance) == True

    def test_bit_pack_dataarray_encoding(self):
        """
        Check DataArray encoding dict doesn't regress.
        """
        data = np.array([-4, 67.12, 80, 1800.2])
        da = xr.DataArray(
            data,
            name="data",
            coords={"x": range(len(data))},
            dims=["x"]
        )
        encoding = bit_pack_dataarray_encoding(da, 16)
        expected = {
            'data': {
                'dtype': 'int16',
                'scale_factor': 0.030758811256559746,
                'add_offset': 898.1,
                '_FillValue': -32768
            }
        }
        assert encoding == expected

    def test_bit_pack_dataarray_encoding_round_trip(self):
        """
        Check DataArray round trip is within tolerance.
        """
        arrays_to_test = [
            np.array([]),
            np.array([np.nan]),
            np.array([np.nan, np.nan, np.nan]),
            np.array([0]),
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
            encoding = bit_pack_dataarray_encoding(da, 16)
            assert self.round_trip_within_tolerance(da, encoding, self.relative_tolerance) == True

    def test_bit_pack_dataarray_encoding_non_finite(self):
        """
        Do not allow bit-packing when data contains infinity.
        """
        arrays_to_test = [
            np.array([np.inf]),
            np.array([0, -np.inf]),
            np.array([0, -np.inf, np.inf]),
        ]
        for data in arrays_to_test:
            da = xr.DataArray(
                data,
                name="data",
                coords={"x": range(len(data))},
                dims=["x"]
            )
            with pytest.raises(AssertionError):
                bit_pack_dataarray_encoding(da, 16)

    def test_bit_pack_dataarray_encoding_degenerate_input(self):
        """
        Test for 0 length input array.
        """
        da = xr.DataArray(
            np.array([]),
            name="data",
            coords={"x": []},
            dims=["x"]
        )
        expected = {
            "data": {
                "dtype": "int16",
                "scale_factor": 1,
                "add_offset": 0,
                "_FillValue": -1
            }
        }
        assert bit_pack_dataarray_encoding(da, 16) == expected

    def test_bit_pack_dataset_encoding(self):
        """
        Test bit-packing on a Dataset of multiple DataArrays.
        """
        np.random.seed(0)
        temperature = 15 + 8 * np.random.randn(2, 2)
        precipitation = 10 * np.random.rand(2, 2)
        lon = [[-99.83, -99.32], [-99.79, -99.23]]
        lat = [[42.25, 42.21], [42.63, 42.59]]

        ds = xr.Dataset(
            data_vars=dict(
                temperature=(["x", "y"], temperature),
                precipitation=(["x", "y"], precipitation),
            ),
            coords=dict(
                lon=(["x", "y"], lon),
                lat=(["x", "y"], lat),
            ),
        )
        expected = {
            'temperature': {
                'dtype': 'int16',
                'scale_factor': 0.0002510535457941546,
                'add_offset': 25.564201630274724,
                '_FillValue': -32768
            },
            'precipitation': {
                'dtype': 'int16',
                'scale_factor': 7.980689171904907e-05,
                'add_offset': 6.577139000604922,
                '_FillValue': -32768
            }
        }

        assert bit_pack_dataset_encoding(ds, 16) == expected