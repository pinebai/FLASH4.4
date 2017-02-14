from StringIO import StringIO

from nose.tools import assert_true, assert_false, assert_equal, assert_raises

from flash.dsl import runtime_parameters


def test_flat_reshape1():
    obs = runtime_parameters._flat_reshape(range(4), (2, 2))
    exp = [(0, 1), (2, 3)]
    assert_equal(obs, exp)


def test_flat_reshape2():
    obs = runtime_parameters._flat_reshape(range(9), (3, 3))
    exp = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    assert_equal(obs, exp)


def test_flat_reshape3():
    obs = runtime_parameters._flat_reshape(range(8), (2, 2, 2))
    exp = [((0, 1), (2, 3)), ((4, 5), (6, 7))]
    assert_equal(obs, exp)


def test_flat_reshape4():
    obs = runtime_parameters._flat_reshape(range(8), (4, 2))
    exp = [(0, 1), (2, 3), (4, 5), (6, 7)]
    assert_equal(obs, exp)


def test_flat_reshape5():
    obs = runtime_parameters._flat_reshape(range(8), (2, 4))
    exp = [(0, 1, 2, 3), (4, 5, 6, 7)]
    assert_equal(obs, exp)


def test_flat_reshape6():
    obs = runtime_parameters._flat_reshape(range(16), (1, 1, 4, 4))
    exp = [(((0, 1, 2, 3), (4, 5, 6, 7), (8, 9, 10, 11), (12, 13, 14, 15)),)]
    assert_equal(obs, exp)


def test_flat_reshape7():
    obs = runtime_parameters._flat_reshape(range(16), (2, 4, 2))
    exp = [((0, 1), (2, 3), (4, 5), (6, 7)), ((8, 9), (10, 11), (12, 13), (14, 15))]
    assert_equal(obs, exp)


def test_nest_flat1():
    a = [((1, 1), 0), ((1, 2), 1), ((2, 1), 2), ((2, 2), 3)]
    obs = runtime_parameters._nest_flat(a)
    exp = [[0, 1], [2, 3]]
    assert_equal(obs, exp)


def test_nest_flat2():
    a = [((1, 1), 0), ((1, 2), 1), ((2, 1), 2), ]
    obs = runtime_parameters._nest_flat(a)
    exp = [[0, 1], [2,]]
    assert_equal(obs, exp)


def test_nest_flat3():
    a = [((1, 1), 0), ((1, 42), 1), ((2, 1), 2), ((2, 2), 3)]
    assert_raises(IndexError, runtime_parameters._nest_flat, a)


def test_nest_flat4():
    a = [((1, 1, 1), 0), 
         ((1, 1, 2), 1), 
         ((1, 2, 1), 2), 
         ((1, 2, 2), 3), 
         ((2, 1, 1), 4), 
         ((2, 1, 2), 5), 
         ((2, 2, 1), 6), 
         ((2, 2, 2), 7), 
         ]
    obs = runtime_parameters._nest_flat(a)
    exp = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
    assert_equal(obs, exp)


def test_nest_flat5():
    a = [((1, 1, 1, 1), 0), 
         ((1, 1, 1, 2), 1), 
         ((1, 1, 1, 3), 2), 
         ((1, 1, 1, 4), 3), 
         ((1, 1, 2, 1), 4), 
         ((1, 1, 2, 2), 5), 
         ((1, 1, 2, 3), 6), 
         ((1, 1, 2, 4), 7), 
         ]
    obs = runtime_parameters._nest_flat(a)
    exp = [[[[0, 1, 2, 3], [4, 5, 6, 7], ]]]
    assert_equal(obs, exp)


def test_nest_flat6():
    a = [((1, 1, 1, 1, 1), 0), 
         ((1, 1, 1, 2, 1), 1), 
         ((1, 1, 1, 3, 1), 2), 
         ((1, 1, 1, 4, 1), 3), 
         ((1, 1, 2, 1, 1), 4), 
         ((1, 1, 2, 2, 1), 5), 
         ((1, 1, 2, 3, 1), 6), 
         ((1, 1, 2, 4, 1), 7), 
         ]
    obs = runtime_parameters._nest_flat(a)
    exp = [[[[[0], [1], [2], [3]], [[4], [5], [6], [7]], ]]]
    assert_equal(obs, exp)


def test_py_form1():
    """Test strings"""
    d = {'a': "'Hello'"}
    obs = runtime_parameters.py_form(d)
    exp = {'a': "Hello"}
    assert_equal(obs, exp)


def test_py_form2():
    """Test strings"""
    d = {'a': '"Hello"'}
    obs = runtime_parameters.py_form(d)
    exp = {'a': "Hello"}
    assert_equal(obs, exp)


def test_py_form3():
    """Test int"""
    d = {'a': "42"}
    obs = runtime_parameters.py_form(d)
    exp = {'a': 42}
    assert_equal(obs, exp)


def test_py_form4():
    """Test floats"""
    d = {'a': "42.0"}
    obs = runtime_parameters.py_form(d)
    exp = {'a': 42.0}
    assert_equal(obs, exp)


def test_py_form5():
    """Test floats"""
    d = {'a': "-4.20e+01"}
    obs = runtime_parameters.py_form(d)
    exp = {'a': -42.0}
    assert_equal(obs, exp)


def test_py_form6():
    """Test bools"""
    d = {'a': ".True."}
    obs = runtime_parameters.py_form(d)
    exp = {'a': True}
    assert_equal(obs, exp)


def test_py_form7():
    """Test bools"""
    d = {'a': ".false."}
    obs = runtime_parameters.py_form(d)
    exp = {'a': False}
    assert_equal(obs, exp)


def test_py_form8():
    """Test indexed"""
    d = {'a_1': "0", 'a_3': "2", 'a_2': "1",}
    obs = runtime_parameters.py_form(d)
    exp = {'a': range(3)}
    assert_equal(obs, exp)


def test_py_form9():
    """Test indexed, not dense."""
    d = {'a_1': "0", 'a_10': "2", 'a_2': "1",}
    obs = runtime_parameters.py_form(d)
    exp = {'a_1': 0, 'a_10': 2, 'a_2': 1,}
    assert_equal(obs, exp)


def test_py_form10():
    """Test indexed"""
    d = {'a_1_2': "0", 'a_2_2': "2", 'a_2_1': "1", 'a_1_1': "10"}
    obs = runtime_parameters.py_form(d)
    exp = {'a': [[10, 0], [1, 2]]}
    assert_equal(obs, exp)


def test_py_form11():
    """Test indexed, not dense."""
    d = {'a_1_2': "0", 'a_2_2': "2", 'a_2_1': "1",}
    obs = runtime_parameters.py_form(d)
    exp = {'a': [[0], [1, 2]]}
    assert_equal(obs, exp)


def test__par_convert_iterable1():
    obs = runtime_parameters._par_convert_iterable('a', [1, 2, 3])
    exp = {'a_1': '1', 'a_2': '2','a_3': '3',}
    assert_equal(obs, exp)
    

def test__par_convert_iterable2():
    obs = runtime_parameters._par_convert_iterable('a', [[1, 2], [3, 4]])
    exp = {'a_1_1': '1', 'a_1_2': '2', 'a_2_1': '3', 'a_2_2': '4',}
    assert_equal(obs, exp)


def test_par_form1():
    obs = runtime_parameters.par_form({'a': [1, 2, 3]})
    exp = {'a_1': '1', 'a_2': '2','a_3': '3',}
    assert_equal(obs, exp)


def test_par_form2():
    obs = runtime_parameters.par_form({'a': 42})
    exp = {'a': '42'}
    assert_equal(obs, exp)


def test_par_form3():
    obs = runtime_parameters.par_form({'a': 'V.v.V'})
    exp = {'a': '"V.v.V"'}
    assert_equal(obs, exp)


def test_par_form4():
    obs = runtime_parameters.par_form({'a': 42.0})
    exp = {'a': "{0:.6e}".format(42.0)}
    assert_equal(obs, exp)


def test_par_form5():
    obs = runtime_parameters.par_form({'a': True})
    exp = {'a': ".true."}
    assert_equal(obs, exp)


def test_par_form6():
    obs = runtime_parameters.par_form({'a': False})
    exp = {'a': ".false."}
    assert_equal(obs, exp)


def test_par_form7():
    obs = runtime_parameters.par_form({'a': 4242424242424242242424242})
    exp = {'a': '4242424242424242242424242'}
    assert_equal(obs, exp)




sample_par = """run_comment = "Simul"
log_file    = "lasslab.log"
basenm      = "lasslab_"

# This particular parfile is used as an example parfile that is
# described in detail in the users guide.
# gr_hypreSolverType = HYPRE_PCG

#########################
#   OUTPUT PARAMETERS   #
#########################

### Checkpoint Options  ###
checkpointFileIntervalTime = 0
checkpointFileIntervalStep = 50

### Plot Options ###
plotFileNumber       = 0
plotFileIntervalStep = 0
plotFileIntervalTime = 1e-9
plot_var_1           = "dens"
plot_var_2           = "depo"
plot_var_3           = "tele"
plot_var_4           = "tion"
plot_var_5           = "ye  "
plot_var_6           = "sumy"
plot_var_6           = "pele"
plot_var_6           = "erad"
plot_var_6           = "targ"

### Restart Options ###
restart              = .false.
checkpointFileNumber = 0

########################################
#                                      #
#     RADIATION/OPACITY PARAMETERS     #
#                                      #
########################################
rt_useMGD       = .true.
rt_mgdNumGroups = 6
rt_mgdBounds_1  = 1.0e-01
rt_mgdBounds_2  = 1.0e+00
rt_mgdBounds_3  = 1.0e+01
rt_mgdBounds_4  = 1.0e+02
rt_mgdBounds_5  = 1.0e+03
rt_mgdBounds_6  = 1.0e+04
rt_mgdBounds_7  = 1.0e+05
rt_mgdFlMode    = "fl_harmonic"
rt_mgdFlCoef    = 1.0

############################
#                          #
#     HYDRO PARAMETERS     #
#                          #
############################

# Use second order hybrid solver with minmod slope limiter. This
# essentially eliminates any Carbuncle instability.

order            = 3        # Interpolation 'order' (first/second/third/fifth order)
slopeLimiter     = "minmod" # Slope limiters (minmod, mc, vanLeer, hybrid, limited)
LimitedSlopeBeta = 1.       # Slope parameter for the "limited" slope by Toro
charLimiting     = .true.   # Characteristic limiting vs. Primitive limiting
use_avisc        = .true.   # use artificial viscosity (originally for PPM)
cvisc            = 0.1      # coefficient for artificial viscosity
use_flattening   = .false.  # use flattening (dissipative) (originally for PPM)
use_steepening   = .false.  # use contact steepening (originally for PPM)
use_upwindTVD    = .false.  # use upwind biased TVD slope for PPM (need nguard=6)
RiemannSolver    = "hll" # Roe, HLL, HLLC, LLF, Marquina, hybrid # initaly at HLL
entropy          = .false.  # Entropy fix for the Roe solver
shockDetect      = .false.  # Shock Detect for '#' numerical stability
gr_hypreMaxIter  = 1000     #  when conjugate-gradient in  HYPRE fails to converg
gr_hypreSolverType = "HYPRE_GMRES" # this seems to be the best solver to use for now

#gr_hypreRelTol = 2.0e-7     #
gr_hypreFloor = 5.0e-7

# Hydro boundary conditions:
xl_boundary_type = "reflect"
xr_boundary_type = "outflow"
yl_boundary_type = "outflow"
yr_boundary_type = "outflow"
zl_boundary_type = "reflect"
zr_boundary_type = "reflect"

### SETUP LASER PULSES ###
ed_numPulses = 1

# Define a luli 2000 1.5ns pulse
ed_numSections_1 = 5
ed_time_1_1  = 0.0e-9
ed_time_1_2  = 0.25e-9
ed_time_1_3  = 1.0e-9
ed_time_1_4  = 1.2e-9
ed_time_1_5  = 1.5e-9

# laser power in Watt
ed_power_1_1  = 0.0
ed_power_1_2  = 4.69724771e+11
ed_power_1_3  = 3.22935780e+11
ed_power_1_4  = 4.69724771e+10
ed_power_1_5  = 0.00"""

def test_load():
    """Integration test for loading par files."""
    f = StringIO(sample_par)
    params = runtime_parameters.load(f)

    assert_equal(params['gr_hypreFloor'], 5.0e-7)
    assert_equal(params['plot_var'][5], "targ")
    assert_equal(params['slopeLimiter'], 'minmod')
    assert_equal(params['order'], 3)
    assert_true(isinstance(params['order'], int))
    assert_false(params['entropy'])
    assert_false(params['shockDetect'])


def test_dump1():
    f = StringIO()
    params = {'b': 10, 'a': True}

    # dump the values
    runtime_parameters.dump(params, f)
    f.pos = 0
    obs = f.read()
    exp = "a = .true.\nb = 10"    

    assert_equal(obs, exp)


def test_load_dump_roundtrip():
    """Integration test for roundtripping par files."""
    f = StringIO(sample_par)
    pydict = runtime_parameters.load(f)

    newf = StringIO()
    runtime_parameters.dump(pydict, newf)
    newf.pos = 0
    obs_dumped = newf.read()

    newf.pos = 0    
    newpydict = runtime_parameters.load(newf)
    assert_equal(set(pydict.keys()), set(newpydict.keys()))

    newnewf = StringIO()
    runtime_parameters.dump(newpydict, newnewf)
    newnewf.pos = 0
    obs_newdumped = newnewf.read()

    assert_equal(obs_dumped, obs_newdumped)
