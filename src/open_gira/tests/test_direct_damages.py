import pytest

from open_gira.direct_damages import ReturnPeriodMap, AqueductFlood


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
