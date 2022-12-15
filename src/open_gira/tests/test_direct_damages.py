import pytest

import numpy as np

from open_gira.direct_damages import (
    ReturnPeriodMap, AqueductFlood, holland_wind_model
)


class TestReturnPeriodMap:
    """
    Test base class rejects instantiation.
    """

    with pytest.raises(TypeError):
        ReturnPeriodMap()


class TestAqueductFlood:
    """
    Test inference of attributes from map name.
    """

    def test_riverine_map_parsing(self):
        name = "inunriver_rcp8p5_00000NorESM1-M_2080_rp00005"
        afm = AqueductFlood(name)
        assert afm.name == name
        assert afm.riverine is True
        assert afm.coastal is False
        assert afm.scenario == "rcp8p5"
        assert afm.model == "00000NorESM1-M"
        with pytest.raises(AttributeError):
            afm.subsidence
        assert afm.year == 2080
        assert afm.return_period_years == 5
        assert afm.annual_probability == 0.2
        assert afm.without_RP == "inunriver_rcp8p5_00000NorESM1-M_2080"

    def test_coastal_map_parsing(self):
        name = "inuncoast_rcp8p5_wtsub_2080_rp1000_0_perc_50"
        afm = AqueductFlood(name)
        assert afm.name == name
        assert afm.riverine is False
        assert afm.coastal is True
        assert afm.scenario == "rcp8p5"
        assert afm.model == "wtsub"
        assert afm.subsidence is True
        assert afm.year == 2080
        assert afm.return_period_years == 1000
        assert afm.annual_probability == 0.001
        assert afm.without_RP == "inuncoast_rcp8p5_wtsub_2080_0_perc_50"

    def test_coastal_map_sea_level_rise_parsing(self):
        assert AqueductFlood("inuncoast_rcp8p5_wtsub_2080_rp1000_0_perc_05").slr_percentile == 5.0
        assert AqueductFlood("inuncoast_rcp8p5_wtsub_2080_rp1000_0_perc_50").slr_percentile == 50.0
        assert AqueductFlood("inuncoast_rcp8p5_wtsub_2080_rp1000_0").slr_percentile == 95.0


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
