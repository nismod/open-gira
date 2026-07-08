"""
Tests for raw CHAZ (netCDF) and Emanuel (MATLAB) ingest, on fabricated
miniature inputs mimicking each source's array layout.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from trackset import physics
from trackset.readers import (
    filter_epoch,
    iter_chaz_ensembles,
    parse_emanuel_arrays,
    process_chaz_ensemble,
    read_chaz_dataset,
)
from trackset.units import MS_PER_KNOT

REFERENCE = pd.Timestamp("1950-01-01")


@pytest.fixture
def chaz_dataset() -> xr.Dataset:
    """
    Two storms, two ensemble members, four timesteps of 6-hourly track data.
    Storm 0 (year 2010, WP) uses all four steps; storm 1 (year 2012, NA) has
    two steps of NaN padding. Ragged arrays are (lifelength, stormID),
    except wind speed which gains a leading ensemble dimension.
    """

    def days_since_reference(timestamp: str, steps: int, pad: int = 0):
        start = (pd.Timestamp(timestamp) - REFERENCE).days
        return [start + 0.25 * step for step in range(steps)] + [np.nan] * pad

    time = np.column_stack(
        [
            days_since_reference("2010-07-01", 4),
            days_since_reference("2012-09-01", 2, pad=2),
        ]
    )
    longitude = np.column_stack(
        [[140.0, 140.5, 141.0, 141.5], [300.0, 300.5] + [np.nan] * 2]
    )
    latitude = np.column_stack([[15.0, 15.5, 16.0, 16.5], [20.0, 20.5] + [np.nan] * 2])
    # knots; second ensemble member is 10kt windier
    wind = np.stack(
        [
            np.column_stack([[40.0, 50.0, 60.0, 55.0], [70.0, 75.0] + [np.nan] * 2]),
            np.column_stack([[50.0, 60.0, 70.0, 65.0], [80.0, 85.0] + [np.nan] * 2]),
        ]
    )

    return xr.Dataset(
        {
            "time": (("lifelength", "stormID"), time),
            "longitude": (("lifelength", "stormID"), longitude),
            "latitude": (("lifelength", "stormID"), latitude),
            "Mwspd": (("ensembleNum", "lifelength", "stormID"), wind),
        },
        coords={
            "lifelength": range(4),
            "stormID": [0, 1],
            "ensembleNum": [0, 1],
        },
    )


class TestChaz:
    def test_iter_ensembles(self, chaz_dataset):
        first, second = iter_chaz_ensembles(chaz_dataset, "CRH", sample=0)

        # NaN padding dropped: 4 + 2 points per ensemble member
        assert len(first) == 6
        assert len(second) == 6

        # track_id encodes genesis method, sample, source year, storm, ensemble
        assert set(first["track_id"]) == {"H_000_2010_00000_00", "H_000_2012_00001_00"}
        assert set(second["track_id"]) == {"H_000_2010_00000_01", "H_000_2012_00001_01"}

        # timestamps reconstructed from days-since-1950
        storm_0 = first[first["storm"] == 0]
        assert storm_0.index[0] == pd.Timestamp("2010-07-01 00:00")
        assert storm_0.index[1] == pd.Timestamp("2010-07-01 06:00")
        assert (storm_0["source_year"] == 2010).all()

        # winds converted from knots
        assert storm_0["max_wind_speed_ms"].iloc[0] == pytest.approx(40.0 * MS_PER_KNOT)

    def test_process_ensemble(self, chaz_dataset):
        (first, _) = iter_chaz_ensembles(chaz_dataset, "CRH", sample=0)
        processed = process_chaz_ensemble(first)

        # basins tagged from position
        by_storm = processed.groupby("storm")["basin_id"].unique()
        assert by_storm[0].tolist() == ["WP"]
        assert by_storm[1].tolist() == ["NA"]

        # radius inferred from wind speed and latitude
        expected_rmw = physics.r_max_willoughby_2004(
            processed["max_wind_speed_ms"], processed["latitude_deg"]
        )
        assert processed["radius_to_max_winds_km"].tolist() == expected_rmw.tolist()

        # inferred pressure is a plausible depression of ambient
        p_env = processed["basin_id"].map(physics.ENV_PRESSURE)
        assert (processed["min_pressure_hpa"] < p_env).all()
        assert (processed["min_pressure_hpa"] > 800.0).all()

    def test_read_dataset_concatenates_ensembles(self, chaz_dataset):
        df = read_chaz_dataset(chaz_dataset, "CRH", sample=0)
        assert len(df) == 12
        assert set(df["ensemble"]) == {0, 1}

    def test_filter_epoch(self):
        df = pd.DataFrame({"source_year": range(2000, 2021)})

        windowed = filter_epoch(df, epoch=2010, half_width_years=5)
        # window bounds are exclusive
        assert windowed["source_year"].min() == 2006
        assert windowed["source_year"].max() == 2014

        # source data must span the whole window
        with pytest.raises(ValueError, match="out of range"):
            filter_epoch(df, epoch=2015, half_width_years=10)


class TestEmanuel:
    @pytest.fixture
    def matlab_arrays(self) -> dict:
        """
        Two tracks, padded to four timesteps: track 0 (year 1999, NA) with
        three observations, track 1 (year 2001, WP) with two. Padding is
        marked by zeros in daystore.
        """

        return {
            "yearstore": np.array([1999, 2001]),
            "monthstore": np.array([[7, 7, 7, 0], [9, 9, 0, 0]]),
            "daystore": np.array([[10, 10, 10, 0], [5, 5, 0, 0]]),
            "hourstore": np.array([[0, 6, 12, 0], [0, 6, 0, 0]]),
            "vstore": np.array([[40.0, 50.0, 45.0, 0.0], [70.0, 75.0, 0.0, 0.0]]),
            "rmstore": np.array([[30.0, 28.0, 29.0, 0.0], [20.0, 18.0, 0.0, 0.0]]),
            "pstore": np.array([[990.0, 985.0, 988.0, 0.0], [960.0, 955.0, 0.0, 0.0]]),
            "longstore": np.array(
                [[300.0, 300.5, 301.0, 0.0], [140.0, 140.5, 0.0, 0.0]]
            ),
            "latstore": np.array([[20.0, 20.5, 21.0, 0.0], [15.0, 15.5, 0.0, 0.0]]),
        }

    def test_parse(self, matlab_arrays):
        df = parse_emanuel_arrays(matlab_arrays)

        # padding dropped: 3 + 2 points
        assert len(df) == 5

        # timestamps assembled from year/month/day/hour arrays
        track_0 = df[df["tc_number"] == 0]
        assert track_0.index[1] == pd.Timestamp("1999-07-10 06:00")
        assert track_0["timestep"].tolist() == [0, 1, 2]

        # winds converted from knots, radii and pressures passed through
        assert track_0["max_wind_speed_ms"].iloc[0] == pytest.approx(40.0 * MS_PER_KNOT)
        assert track_0["radius_to_max_winds_km"].tolist() == [30.0, 28.0, 29.0]
        assert track_0["min_pressure_hpa"].tolist() == [990.0, 985.0, 988.0]

        # basin re-tagged from position, track_id built from it
        assert set(df["basin_id"]) == {"NA", "WP"}
        assert set(df["track_id"]) == {"NA_1999_0000", "WP_2001_0001"}

        assert (df["sample"] == 0).all()
