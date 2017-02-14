from StringIO import StringIO

from nose.tools import assert_true, assert_false, assert_equal, assert_raises
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from flash import output


sample_laser_dat = """#Step     Time (s)     dt (s)      Energy in (erg)    Energy out (erg)       dE in (erg)        dE out (erg)
#
    1   1.0000E-13  1.0000E-13    4.5000000000E+05    4.4196714106E+05    4.5000000000E+05    4.4196714106E+05
    2   1.6185E-13  6.1850E-14    7.2832345640E+05    7.2028436116E+05    2.7832345640E+05    2.7831722010E+05
    3   1.6686E-13  5.0150E-15    7.5089111325E+05    7.4285136089E+05    2.2567656850E+04    2.2566999721E+04
    4   1.7284E-13  5.9763E-15    7.7778444062E+05    7.6974356297E+05    2.6893327377E+04    2.6892202087E+04
    5   1.7933E-13  6.4887E-15    8.0698368656E+05    7.9894088309E+05    2.9199245933E+04    2.9197320114E+04
    6   1.8541E-13  6.0787E-15    8.3433803222E+05    8.2629223557E+05    2.7354345660E+04    2.7351352485E+04
"""

def test_load_laser_dat():
    laser_dat = StringIO(sample_laser_dat)
    obs = output.load_laser_dat(laser_dat)
    raw = [tuple(map(float, sld.split())) for sld in sample_laser_dat.split('\n')[2:-1]]
    exp = np.array(raw, dtype=[('step', np.int32),
                               ('t', np.float64),
                               ('dt', np.float64),
                               ('E_in', np.float64),
                               ('E_out', np.float64),
                               ('dE_in', np.float64),
                               ('dE_out', np.float64),
                               ])
    assert_array_equal(obs, exp)
