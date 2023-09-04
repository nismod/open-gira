import numpy as np

from open_gira.geodesic import bearing_and_great_circle_distance


class TestBearingAndGreatCircleDistance:
    """
    Test function which computes forward bearing and great-circle distance
    (meters) between points.
    """

    def test_point_pair(self):
        lon0 = -81.78
        lat0 = 24.56
        lon1 = -81.09
        lat1 = 24.71
        pts = [lon0, lat0, lon1, lat1]

        assert np.allclose(
            bearing_and_great_circle_distance(*pts),
            (76.40617746089684, 71707.82335842062)
        )

    def test_coincident_points(self):
        lon0 = 0.9764181
        lat0 = 50.9232914
        lon1 = 0.9764181
        lat1 = 50.9232914
        pts = [lon0, lat0, lon1, lat1]

        # bearing for coincident points is ill-defined, do not test
        _, distance = bearing_and_great_circle_distance(*pts)
        assert distance == 0

    def test_antipodes(self):
        lon0 = 10
        lat0 = 10
        lon1 = -170
        lat1 = -10
        pts = [lon0, lat0, lon1, lat1]

        # bearing for antipodes is ill-defined, do not test
        _, distance = bearing_and_great_circle_distance(*pts)

        # earth circumference ~40M meters
        assert np.allclose(distance, 20015086.606149975)

    def test_point_arrays(self):
        # three random points
        pts = [
            np.array([-78.85524587, 45.46099115, -28.3258711]),
            np.array([28.92715285, 79.39721556, -44.00797675]),
            np.array([-156.00968598, -84.57826798, 158.38471011]),
            np.array([-54.79011992, 15.13498674, 72.21865346])
        ]

        assert np.allclose(
            bearing_and_great_circle_distance(*pts),
            (
                np.array([-144.11891165, -48.30133907, -4.30392061]),
                np.array([11835578.78885653, 9097391.66244342, 16857981.86371326])
            )
        )
