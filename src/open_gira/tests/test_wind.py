import numpy as np

from open_gira.wind import holland_wind_model, advective_vector, sigmoid_decay, estimate_wind_field


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


class TestSigmoidDecay:
    """Test sigmoid/logistic decay function."""

    def test_sigmoid_decay(self):
        expected = [
            0.98201379, 0.97340301, 0.96083428, 0.94267582, 0.9168273 ,
            0.88079708, 0.83201839, 0.76852478, 0.68997448, 0.59868766,
            0.5       , 0.40131234, 0.31002552, 0.23147522, 0.16798161,
            0.11920292, 0.0831727 , 0.05732418, 0.03916572, 0.02659699,
            0.01798621
        ]

        np.testing.assert_allclose(
            expected,
            sigmoid_decay(np.arange(0, 1001, 50), 500, 0.004),
            rtol=1E-6
        )


class TestEstimateWindField:
    """
    Test function which computes fields of advective and rotational vectors and combines them.

    Testing a transect in longitude either side of storm eye, with single valued latitude.
    """

    def test_estimate_wind_field_rotation_only(self):
        # stationary, rotating with maximum velocity of 80ms-1 at 10,000m radius
        # should rotate anticlockwise in northern hemisphere
        result = estimate_wind_field(np.linspace(-30.5, -29.5, 9), 15, -30, 15, 10000, 80, 92000, 100000, 0, 0)
        np.testing.assert_allclose(
            np.abs(result),
            np.array(
                [[16.58826006, 23.81669229, 38.22875421, 72.35172612, np.nan,
                72.35172612, 38.22875421, 23.81669229, 16.58826006]]
            )
        )
        np.testing.assert_allclose(
            np.angle(result, deg=True),
            np.array(
                # degrees CCW from positive real axis (mathematical convention)
                [[-89.93529486, -89.95147127, -89.96764757, -89.9838238 ,
                np.nan,  89.9838238 ,  89.96764757,  89.95147127, 89.93529486]]
            )
        )

    def test_estimate_wind_field_rotation_only_southern_hemisphere(self):
        # stationary, rotating with maximum velocity of 80ms-1 at 10,000m radius
        # should rotate clockwise in southern hemisphere
        result = estimate_wind_field(np.linspace(-30.5, -29.5, 9), -15, -30, -15, 10000, 80, 92000, 100000, 0, 0)
        np.testing.assert_allclose(
            np.abs(result),
            np.array(
                [[16.58826006, 23.81669229, 38.22875421, 72.35172612, np.nan,
                72.35172612, 38.22875421, 23.81669229, 16.58826006]]
            )
        )
        np.testing.assert_allclose(
            np.angle(result, deg=True),
            np.array(
                # flipped sign w.r.t. test_estimate_wind_field_rotation_only
                # degrees anticlockwise from positive real axis (mathematical convention)
                [[ 89.93529486,  89.95147127,  89.96764757,  89.9838238 , np.nan,
                -89.9838238 , -89.96764757, -89.95147127, -89.93529486]]
            )
        )

    def test_estimate_wind_field_advection_only(self):
        # storm moving N (heading 0 degrees) at 10ms-1, no pressure drop and no rotation (what a perverse storm)
        result = estimate_wind_field(np.linspace(-30.5, -29.5, 9), 15, -30, 15, 10000, 1E-6, 99999, 100000, 0, 10)
        # recovering alpha parameter, ~0.56 factor reduction from eye speed
        np.testing.assert_allclose(
            np.abs(result) / 10,
            np.array(
                [[0.54467011, 0.5461928 , 0.5475677 , 0.54880849, np.nan,
                0.54880849, 0.5475677 , 0.5461928 , 0.54467011]]
            )
        )
        # recovering beta parameter, 19.2 degrees over-rotation
        np.testing.assert_allclose(
            np.angle(result, deg=True) - 90,
            np.array(
                # degrees CCW from positive real axis (mathematical convention)
                [[19.2, 19.2, 19.2, 19.2, np.nan, 19.2, 19.2, 19.2, 19.2]]
            )
        )

    def test_estimate_wind_field_rotation_and_advection(self):
        # storm moving N at 10ms-1, rotating at 80ms-1 with RMW=10km, 100mbar pressure drop
        result = estimate_wind_field(np.linspace(-30.5, -29.5, 9), 15, -30, 15, 10000, 80, 90000, 100000, 0, 10)
        # greater speeds on east of eye, where advection and rotation combine
        np.testing.assert_allclose(
            np.abs(result),
            np.array(
                [[24.45427248, 31.92847355, 44.33622784, 65.62543488, np.nan,
                75.9877587 , 54.67172047, 42.23278268, 34.72300576]]
            )
        )
        np.testing.assert_allclose(
            np.angle(result, deg=True),
            np.array(
                # degrees CCW from positive real axis (mathematical convention)
                [[-94.12223634, -93.16869142, -92.29164582, -91.55850813,
                np.nan,  91.34593479,  91.8582487 ,  92.39504393, 92.90189219]]
            )
        )