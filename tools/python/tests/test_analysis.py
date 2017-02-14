from StringIO import StringIO

from nose.tools import assert_true, assert_false, assert_equal, assert_raises
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from flash import analysis


def test_shock_detect1():
    r = np.arange(10.0)
    v = np.zeros(10, dtype=float)
    v[3] = 1.0

    obs_r, obs_v, obs_i = analysis.shock_detect(r, v)
    exp_r, exp_v, exp_i = 3.0, 1.0, 3

    assert_equal(obs_r, exp_r)
    assert_equal(obs_v, exp_v)
    assert_equal(obs_i, exp_i)


def test_shock_detect2():
    r = np.arange(10.0)
    v = np.zeros(10, dtype=float)
    assert_raises(ValueError, analysis.shock_detect, r, v)
