import numpy as np
import pytest

from trackset import physics


def test_r_max_willoughby_shrinks_with_wind_speed():
    weak = physics.r_max_willoughby_2004(20.0, 20.0)
    strong = physics.r_max_willoughby_2004(70.0, 20.0)
    assert strong < weak
    # spot value: 46.29 * exp(-0.0153 * 20 + 0.0166 * 20)
    assert weak == pytest.approx(46.29 * np.exp(0.026))


def test_r_max_willoughby_hemisphere_symmetry():
    north = physics.r_max_willoughby_2004(np.array([30.0]), np.array([15.0]))
    south = physics.r_max_willoughby_2004(np.array([30.0]), np.array([-15.0]))
    assert north == pytest.approx(south)


def test_coriolis():
    assert physics.coriolis(0.0) == pytest.approx(0.0)
    # |2 * omega * sin(90)| ~ 1.4544e-4 rad/s at the pole
    assert physics.coriolis(90.0) == pytest.approx(1.4544e-4, rel=1e-3)
    assert physics.coriolis(-90.0) == pytest.approx(physics.coriolis(90.0))


def test_b_parameter_in_plausible_range():
    r_max_m = physics.r_max_willoughby_2004(50.0, 20.0) * 1_000
    b = physics.b_vickery_wadhera_2008(20.0, r_max_m)
    # Holland's B is typically ~1-2.5
    assert 0.5 < b < 2.5


class TestPMinHolland:
    def test_deeper_low_for_stronger_storm(self):
        p_env = physics.ENV_PRESSURE["NA"]
        args = dict(phi=20.0, r_max=30_000.0)
        weak = physics.p_min_holland_1980(p_env, v_max=25.0, **args)
        strong = physics.p_min_holland_1980(p_env, v_max=60.0, **args)
        assert strong < weak < p_env

    def test_implausible_pressures_are_nan(self):
        # absurd wind speed drives estimate below 800 hPa
        p_min = physics.p_min_holland_1980(
            1010.0, v_max=200.0, r_max=30_000.0, phi=20.0
        )
        assert np.isnan(p_min)
