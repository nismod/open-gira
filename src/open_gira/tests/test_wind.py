import numpy as np

from open_gira.wind import holland_wind_model, advective_vector, rotational_field


class TestHollandWindModel:
    """Test function which computes the radial profile of storm wind speed"""

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


class TestAdvectiveVector:
    """Test function which computes advective vector from storm velocity"""

    def test_simple_model(self):
        # alpha=1, beta=0, i.e. advective component _is_ eye velocity
        np.testing.assert_allclose(
            advective_vector(0, 10, 1, 1, 0),
            0 + 10j
        )

    def test_alpha_scaling(self):
        np.testing.assert_allclose(
            advective_vector(0, 10, 1, 0.5, 0),
            0 + 5j
        )

    def test_beta_rotation_northern_hemisphere(self):
        # northern hemisphere, storms rotate anticlockwise
        # therefore, beta=90 should rotate vector from N to W
        np.testing.assert_allclose(
            advective_vector(0, 10, 1, 1, 90),
            -10 + 0j
        )

    def test_beta_rotation_southern_hemisphere(self):
        # southern hemisphere, storms rotate clockwise
        # therefore, beta=90 should rotate vector from N to E
        np.testing.assert_allclose(
            advective_vector(0, 10, -1, 1, 90),
            10 + 0j
        )

    def test_default_argument_regression(self):
        # ensure alpha and beta aren't changed by accident
        np.testing.assert_allclose(
            advective_vector(0, 10, 1),
            -1.84165322 + 5.28850767j
        )


class TestRotationalField:
    """
    Test function which computes field of rotational vectors, that is intensity
    from the Holland model, and rotation direction from hemisphere.
    """

    def test_rotational_field(self):
        eye_lon = -30
        lon_width = 5
        max_lon = eye_lon + (lon_width - 1) / 2
        min_lon = eye_lon - (lon_width - 1) / 2
        lon_arr = np.linspace(min_lon, max_lon, lon_width)
        eye_lat = 15

        expected = np.array(
            [[0.00063 -0.139398j, 0.029922-13.247382j, np.nan + np.nan * 1j, 0.029922+13.247382j, 0.00063  +0.139398j]]
        )
        # test for regression
        north = rotational_field(lon_arr, np.array([eye_lat]), eye_lon, eye_lat, 50_000, 100, 95_000, 100_000)
        np.testing.assert_allclose(
            north,
            expected,
            rtol=1E-5
        )

        # check that the direction of rotation is reversed in the southern hemisphere
        south = rotational_field(lon_arr, np.array([-eye_lat]), eye_lon, -eye_lat, 50_000, 100, 95_000, 100_000)

        np.testing.assert_allclose(
            np.angle(north),
            -1 * np.angle(south)
        )
