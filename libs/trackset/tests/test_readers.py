import textwrap

import geopandas as gpd
import pandas as pd
import pytest

from trackset import (
    read_chaz,
    read_emanuel,
    read_ibtracs,
    read_iris,
    read_storm,
)
from trackset.units import MS_PER_KNOT, TEN_MINUTE_TO_ONE_MINUTE_WIND_FACTOR


class TestReadStorm:
    # year,month,tc_number,timestep,basin_id,lat,lon,min_pressure_hpa,
    # max_wind_speed_ms,radius_to_max_winds_km,category,landfall,distance_to_land_km
    # includes a duplicate timestep (year 0, tc 1, timestep 1) and a longitude
    # beyond 180
    CSV = textwrap.dedent(
        """
        0,1,1,0,4,-18.0,178.0,995.0,22.0,30.0,0,0,100.0
        0,1,1,1,4,-18.1,179.2,990.0,25.0,28.0,0,0,90.0
        0,1,1,1,4,-18.1,179.2,990.0,25.0,28.0,0,0,90.0
        0,1,1,2,4,-18.2,182.4,985.0,28.0,26.0,0,0,80.0
        1,2,1,0,4,-15.0,175.0,1000.0,18.0,40.0,0,0,200.0
        1,2,1,1,4,-15.5,176.0,998.0,20.0,38.0,0,0,150.0
        """
    ).strip()

    @pytest.fixture
    def storm_file(self, tmp_path):
        path = tmp_path / "STORM_DATA_IBTRACS_SP_1000_YEARS_2.csv"
        path.write_text(self.CSV)
        return path

    def test_read(self, storm_file):
        ts = read_storm([storm_file], sample=2)

        assert ts.synthetic_time
        assert ts.years == 1000.0
        assert ts.n_tracks == 2

        # duplicate timestep dropped
        assert len(ts.data) == 5

        # integer basin decoded, sample recorded, year offset by 1000 * sample
        assert set(ts.data["basin_id"]) == {"SP"}
        assert set(ts.data["sample"]) == {2}
        assert set(ts.data["year"]) == {2000, 2001}
        assert set(ts.data["track_id"]) == {"SP_2_2000_1", "SP_2_2001_1"}

        # longitudes wrapped into [-180, 180)
        assert ts.data.geometry.x.max() < 180.0
        assert ts.data.geometry.x.min() == pytest.approx(-177.6)

        # winds rescaled from 10-minute to 1-minute sustained
        assert ts.data["max_wind_speed_ms"].iloc[0] == pytest.approx(
            22.0 / TEN_MINUTE_TO_ONE_MINUTE_WIND_FACTOR
        )


class TestReadIbtracs:
    CSV = textwrap.dedent(
        """
        SID,SEASON,NUMBER,BASIN,NAME,ISO_TIME,LAT,LON,DIST2LAND,USA_WIND,USA_PRES,USA_RMW,TOKYO_WIND
        ,Year,,,,,degrees_north,degrees_east,km,kts,mb,nmile,kts
        2017260N12310,2017,25,NA,MARIA,2017-09-19 00:00:00,15.3,-61.1,20,50,985,20,53.3
        2017260N12310,2017,25,NA,MARIA,2017-09-19 03:00:00,15.7,-61.6,0,55,980,18,
        2017260N12310,2017,25,NA,MARIA,2017-09-19 06:00:00,16.1,-62.1,10,60,975, ,
        """
    ).strip()

    @pytest.fixture
    def ibtracs_file(self, tmp_path):
        path = tmp_path / "ibtracs.since1980.list.v04r00.csv"
        path.write_text(self.CSV)
        return path

    def test_read(self, ibtracs_file):
        ts = read_ibtracs(ibtracs_file)

        assert ts.source == "IBTrACS"
        assert not ts.synthetic_time
        assert ts.years == 1.0  # single season observed
        assert ts.n_tracks == 1

        # third point has no radius report from any agency: dropped
        assert len(ts.data) == 2
        assert ts.data["timestep"].tolist() == [0, 1]

        # 'NA' basin must survive as the North Atlantic, not become NaN
        assert set(ts.data["basin_id"]) == {"NA"}
        assert set(ts.data["name"]) == {"MARIA"}

        # first point: TOKYO 53.3 kt normalises to (53.3 - 23.3) / 0.6 = 50 kt,
        # averaging with USA's 50 kt gives 50 kt exactly
        assert ts.data["max_wind_speed_ms"].iloc[0] == pytest.approx(50.0 * MS_PER_KNOT)

        # landfall where distance to land is zero
        assert ts.data["landfall"].tolist() == [False, True]

        # real timestamps preserved
        assert ts.data.index[0] == pd.Timestamp("2017-09-19 00:00:00")


class TestReadIris:
    TXT = textwrap.dedent(
        """
        IRIS synthetic tracks
        tcid year tc month timestep lon lat vmax pmin rmw r18
        WP0001 5 1 7 0 140.0 15.0 25.0 990.0 30.0 100.0
        WP0001 5 1 7 1 141.0 15.5 30.0 985.0 28.0 110.0
        WP0002 6 1 8 0 210.0 20.0 40.0 970.0 25.0 120.0
        WP0002 6 1 8 1 211.0 20.5 45.0 965.0 22.0 130.0
        """
    ).strip()

    @pytest.fixture
    def iris_file(self, tmp_path):
        path = tmp_path / "IRIS_WP_1000Y_n3.txt"
        path.write_text(self.TXT)
        return path

    def test_read(self, iris_file):
        ts = read_iris([iris_file])

        assert ts.synthetic_time
        assert ts.years == 1000.0
        assert ts.n_tracks == 2

        # sample and basin parsed from the file name; year offset by sample
        assert set(ts.data["sample"]) == {3}
        assert set(ts.data["basin_id"]) == {"WP"}
        assert set(ts.data["year"]) == {3005, 3006}
        assert set(ts.data["track_id"]) == {"WP_3_3005_1", "WP_3_3006_1"}

        # longitudes wrapped
        assert ts.data.geometry.x.tolist() == [140.0, 141.0, -150.0, -149.0]


@pytest.fixture
def preparsed_file(tmp_path):
    """A tiny pre-parsed (CHAZ/Emanuel-style) tabular track GeoParquet."""

    df = pd.DataFrame(
        {
            "sample": [0, 0, 0, 0],
            "year": [2000, 2000, 2001, 2001],
            "tc_number": [1, 1, 1, 1],
            "timestep": [0, 1, 0, 1],
            "max_wind_speed_ms": [25.0, 30.0, 40.0, 45.0],
            "min_pressure_hpa": [990.0, 985.0, 970.0, 965.0],
            "radius_to_max_winds_km": [30.0, 28.0, 25.0, 22.0],
        },
        index=pd.date_range("2000-01-01", periods=4, freq="6h"),
    )
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy([220.0, 221.0, 150.0, 151.0], [10, 11, 12, 13]),
        crs=4326,
    )
    path = tmp_path / "raw_tracks.gpq"
    gdf.to_parquet(path)
    return path


def test_read_chaz(preparsed_file):
    ts = read_chaz(
        preparsed_file, source="CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050", years=500.0
    )
    assert ts.years == 500.0
    assert set(ts.data["track_id"]) == {"S000Y2000N001", "S000Y2001N001"}
    # longitudes wrapped
    assert ts.data.geometry.x.tolist() == [-140.0, -139.0, 150.0, 151.0]


def test_read_emanuel(preparsed_file):
    ts = read_emanuel(preparsed_file, source="emanuel_ssp-585", years=200.0)
    assert ts.n_tracks == 2
    assert ts.annual_frequency == pytest.approx(0.01)
