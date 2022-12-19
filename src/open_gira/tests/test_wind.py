import numpy as np

from open_gira.wind import holland_wind_model


class TestHollandWindModel:

    def test_delta_P_zero(self):
        """
        No pressure drop, not really a storm...
        With no pressure drop, we have a divide by zero and therefore NaN wind speed
        """
        args = [1_000, 50, 100000, 100000, np.array([10_000]), 30]
        expected_result = np.array([np.nan])
        np.testing.assert_allclose(holland_wind_model(*args), expected_result)

    def test_radius_zero(self):
        """
        Winds at the centre of the eye
        Radius zero results in a divide by zero and therefore NaN wind speed
        """
        args = [1_000, 50, 98000, 101000, np.array([0]), 30]
        expected_result = np.array([np.nan])
        np.testing.assert_allclose(holland_wind_model(*args), expected_result)

    def test_1D(self):
        """1D array of distances to calculate wind speeds for"""
        args = [1_000, 50, 98000, 101000, np.linspace(10, 10_000, 4), 30]
        expected_result = np.array(
            [ 0.      , 16.289285,  6.453159,  3.577267]
        )
        np.testing.assert_allclose(holland_wind_model(*args), expected_result, rtol=1E-6)

    def test_2D(self):
        """2D array of distances to calculate wind speeds for"""
        X, Y = np.meshgrid(np.linspace(10, 10_000, 3), np.linspace(10, 10_000, 3))
        args = [1_000, 50, 98000, 101000, np.sqrt(X**2 + Y**2), 30]
        expected_result = np.array([
            [0.      , 9.564629, 3.577264],
            [9.564629, 5.936819, 3.002372],
            [3.577264, 3.002372, 2.01601 ]
        ])
        np.testing.assert_allclose(holland_wind_model(*args), expected_result, rtol=1E-6)
