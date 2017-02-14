"""Helper functions for reading and modifying flash output files. For more 
infromation on using this module to make figures with matplotlib and yt, please
refer to the following presentation:
http://flash.uchicago.edu/site/flashcode/user_support/yt_flash_informal.pdf
"""

import collections

import numpy as np
from scipy.interpolate import griddata

try:
    import yt.mods
except ImportError:
    pass

from . import analysis
from . import yt_derived_fields


field_labels = {
    'dens': 'Density [g/cm$^3$]',
    'pres': 'Pressure [J/cm$^3$]',
    'edens': 'Electron Density [cm$^{-3}$]',
    'nele': 'Electron Density [cm$^{-3}$]',
    'nion': 'Ion Density [cm$^{-3}$]',
    'targ': 'Target Material Fraction',
    'tele': 'Electron Temperature ($T_e$) [K]',
    'tion': 'Ion Temperature ($T_i$) [K]',
    'pele': 'Electron Pressure ($P_e$) [J/cm$^3$]',
    'trad': 'Radiation Temperarture [K]',
    'prad': 'Radiation Pressure [J/cm$^3$]',
    'temp': 'Temperature ($T$) [K]',
    'depo': 'Energy Deposition [ergs/g]',
    'ye': '$Y_e$',
    'abar': r'$\bar{A}$',
    'velx': r'Velocity $\hat{x}$ [cm/s]',
    'vely': r'Velocity $\hat{y}$ [cm/s]',
    'velz': r'Velocity $\hat{z}$ [cm/s]',
    'velo': r'Flow Velocity [cm/s]',
    'velr_rz': r'Velocity $\hat{r}$ [cm/s]',
    'velz_rz': r'Velocity $\hat{z}$ [cm/s]',
    'velt_rz': r'Velocity $\hat{\theta}$ [cm/s]',
    }

field_scale_names = {
    'dens': 'log10',
    'pres': 'log10',
    'edens': 'log10',
    'nele': 'log10',
    'nion': 'log10',
    'targ': 'linear',
    'tele': 'log10',
    'pele': 'log10',
    'trad': 'log10',
    'prad': 'log10',
    'temp': 'log10',
    'depo': 'log10',
    'ye': 'linear',
    'abar': 'linear',
    'velx': 'linear',
    'vely': 'linear',
    'velz': 'linear',
    'velo': 'linear',
    'velr_rz': 'linear',
    'velz_rz': 'linear',
    'velt_rz': 'linear',
    }

_linear = lambda *args: args[0] if 1 == len(args) else args

def _ax_decade_ticks(ax, xory):
    get_ticks = ax.get_xticks if xory == 'x' else ax.get_yticks
    set_ticks = ax.set_xticks if xory == 'x' else ax.set_yticks
    set_labels = ax.set_xticklabels if xory == 'x' else ax.set_yticklabels
    ticks = get_ticks()
    ord_min = int(np.trunc(ticks[0]))
    ord_max = int(np.trunc(ticks[-1]))
    decades = np.arange(ord_min, ord_max + 1)
    ticklabels = ["$10^{{{0}}}$".format(d) for d in decades]
    set_ticks(decades)
    set_labels(ticklabels)
    return ax

def _cb_decade_ticks(cb):
    ord_min = int(np.trunc(cb.vmin))
    ord_max = int(np.trunc(cb.vmax))
    decades = np.arange(ord_min, ord_max + 1)
    ticklabels = ["$10^{{{0}}}$".format(d) for d in decades]
    cb.set_ticks(decades)
    cb.set_ticklabels(ticklabels)
    return cb

field_axis = {}
field_scales = {}
field_colorbars = {}
for _field, _scale_name in field_scale_names.items():
    if _scale_name == 'log10':
        field_axis[_field] = _ax_decade_ticks
        field_scales[_field] = np.log10 
        field_colorbars[_field] = _cb_decade_ticks
    else:
        field_axis[_field] = _linear
        field_scales[_field] = _linear
        field_colorbars[_field] = _linear
del _field, _scale_name

def _to_hashable(d):
    """Returns a version of the dictionary whose items are hashable."""
    hd = {}
    for key in d:
        if isinstance(d[key], collections.MutableSequence):
            hd[key] = tuple(d[key])
        elif isinstance(d[key], collections.MutableSet):
            hd[key] = frozenset(d[key])
        else:
            hd[key] = d[key]
    return hd
    


slice_cache = {}

def slice(axis, coord, field, pf, bounds=None, resolution=600, method="nearest", **kwargs):
    """Grabs the slice of a certain field (parameter, data) from the flash
    output file pf.  This slice is performed at the coord along the axis.

    This function requires both scipy and yt and returns x, y, and z data 
    that is suitable for plotting with matplotlib.imshow().
    
    Parameters
    ----------
    axis : int
        The axis along which to slice.  Can be 0, 1, or 2 for x, y, z.
    coord : float
        The coordinate along the axis at which to slice.  This is in
        "domain" coordinates.
    field : str
        A field to retrieve, e.g. 'dens'.
    pf : str or yt.data_objects.static_output
        The checkpoint or plot file, may be either a path to the file 
        or the yt object returned by pf.mods.load().
    bounds : len-4 sequence of floats or None, optional
        This defines the area within the domain that should be sliced.
        By default, the whole domain is taken.  Applies the bounds in
        axes order, ie if axis=2 then bounds=[xmin, xmax, ymin, ymax].
    resolution : int or len-2 sequence of ints, optional
        If a list of intergers, this is the number of points (pixels) 
        along each dimension.  If a single integer, then this is the 
        resolution the longest axis.  The resolution on the other axis
        is then calculated to preserve the aspect ratio.
    method : str, optional
        Interpolation method flag, passed directly down to the function 
        scipy.interpolate.griddata(). 
    kwargs : optional
        All other keyword arguments are passes to the slice() method
        on the pf.h object from yt.

    Returns
    -------
    xdat : 2d numpy array
        Meshed data for x.
    ydat : 2d numpy array
        Meshed data for y.
    zdat : 2d numpy array
        Meshed data for z.

    """
    # open the file if we got a path.
    opened_here = False
    if isinstance(pf, basestring):
        pf = yt.mods.load(pf)
        opened_here = True

    # make kwargs hashable
    kwargs = _to_hashable(kwargs)

    # get the slice object
    slice_key = (pf, axis, coord) + tuple(sorted(kwargs.items()))
    if slice_key in slice_cache:
        amr_slice = slice_cache[slice_key]
    else:
        amr_slice = pf.h.slice(axis, coord, field, **kwargs)
        slice_cache[slice_key] = amr_slice

    # get points
    x = amr_slice['x']
    y = amr_slice['y']
    z = amr_slice['z']
    points = [(y, z), (x, z), (x, y)][axis]

    # get full domain
    if bounds is None:
        bounds = (points[0][0], points[0][-1], points[1][0], points[1][-1])

    # get proper resolution
    if isinstance(resolution, int):
        lengths = (bounds[1] - bounds[0], bounds[3] - bounds[2])
        min_ind, min_len = min(enumerate(lengths), key=lambda x: x[1])
        resolution = [resolution, resolution]
        resolution[min_ind] = int(resolution[min_ind] * min_len / lengths[(min_ind + 1) % 2])
    elif isinstance(resolution, collections.Sequence) and 2 == len(resolution):
        pass
    else: 
        raise ValueError("resolution must be int or length 2 sequence of ints.")

    # get meshed data
    aa = np.linspace(bounds[0], bounds[1], resolution[0])
    bb = np.linspace(bounds[2], bounds[3], resolution[1])
    adat, bdat = np.meshgrid(aa, bb)
    cdat = griddata(points, amr_slice[field], (adat, bdat), method=method)

    # close the file if we opened it
    if opened_here:
        del pf

    # return in the proper order
    dat = [(cdat, adat, bdat), (adat, cdat, bdat), (adat, bdat, cdat)]
    return dat[axis]


def slice_gradient(axis, coord, field, pf, bounds=None, resolution=600, method="nearest", **kwargs):
    """Grabs the gradient of a slice of a certain field (parameter, data) 
    from the flash output file pf.  This slice is performed at the coord 
    along the axis.

    This function requires both scipy and yt and returns x, y, and z data 
    that is suitable for plotting with matplotlib.imshow().
    
    Parameters
    ----------
    axis : int
        The axis along which to slice.  Can be 0, 1, or 2 for x, y, z.
    coord : float
        The coordinate along the axis at which to slice.  This is in
        "domain" coordinates.
    field : str
        A field to retrieve, e.g. 'dens'.
    pf : str or yt.data_objects.static_output
        The checkpoint or plot file, may be either a path to the file 
        or the yt object returned by pf.mods.load().
    bounds : len-4 sequence of floats or None, optional
        This defines the area within the domain that should be sliced.
        By default, the whole domain is taken.  Applies the bounds in
        axes order, ie if axis=2 then bounds=[xmin, xmax, ymin, ymax].
    resolution : int or len-2 sequence of ints, optional
        If a list of intergers, this is the number of points (pixels) 
        along each dimension.  If a single integer, then this is the 
        resolution the longest axis.  The resolution on the other axis
        is then calculated to preserve the aspect ratio.
    method : str, optional
        Interpolation method flag, passed directly down to the function 
        scipy.interpolate.griddata(). 
    kwargs : optional
        All other keyword arguments are passes to the slice() method
        on the pf.h object from yt.

    Returns
    -------
    adat : 2d numpy array
        Meshed data for the first non-axis axis, eg x when axis=z.
    bdat : 2d numpy array
        Meshed data for the second non-axis axis, eg y when axis=z.
    dfdadat : 2d numpy array
        Gradient of field along a.
    dfdbdat : 2d numpy array
        Gradient of field along b.
    magdat : 2d numpy array
        Magnitude of the gradient of field.

    """
    # open the file if we got a path.
    opened_here = False
    if isinstance(pf, basestring):
        pf = yt.mods.load(pf)
        opened_here = True

    # make kwargs hashable
    kwargs = _to_hashable(kwargs)

    # get the slice object
    slice_key = (pf, axis, coord) + tuple(sorted(kwargs.items()))
    if slice_key in slice_cache:
        amr_slice = slice_cache[slice_key]
    else:
        amr_slice = pf.h.slice(axis, coord, field, **kwargs)
        slice_cache[slice_key] = amr_slice

    # get a, b, f
    x = amr_slice['x']
    y = amr_slice['y']
    z = amr_slice['z']
    a, b = [(y, z), (x, z), (x, y)][axis]

    dx = amr_slice['dx']
    dy = amr_slice['dy']
    dz = amr_slice['dz']
    da, db = [(dy, dz), (dx, dz), (dx, dy)][axis]

    f = amr_slice[field]
    amr_len = len(f)

    # calc gradient along a
    dfda = np.empty_like(f)
    bunique = np.unique(b)
    for bu in bunique:
        mask = bu == b
        dfda[mask] = np.gradient(f[mask]) / da[mask]

    # calc gradient along b
    dfdb = np.empty_like(f)
    aunique = np.unique(a)
    for au in aunique:
        mask = au == a
        dfdb[mask] = np.gradient(f[mask]) / db[mask]

    # calc gradient magnitude
    mag = np.sqrt(dfda**2 + dfdb**2)

    # get full domain
    if bounds is None:
        bounds = (a[0], a[-1], b[0], b[-1])

    # get proper resolution
    if isinstance(resolution, int):
        lengths = (bounds[1] - bounds[0], bounds[3] - bounds[2])
        min_ind, min_len = min(enumerate(lengths), key=lambda x: x[1])
        resolution = [resolution, resolution]
        resolution[min_ind] = int(resolution[min_ind] * min_len / lengths[(min_ind + 1) % 2])
    elif isinstance(resolution, collections.Sequence) and 2 == len(resolution):
        pass
    else: 
        raise ValueError("resolution must be int or length 2 sequence of ints.")

    # get meshed data
    aa = np.linspace(bounds[0], bounds[1], resolution[0])
    bb = np.linspace(bounds[2], bounds[3], resolution[1])
    adat, bdat = np.meshgrid(aa, bb)

    dfdadat = griddata((a, b), dfda, (adat, bdat), method=method)
    dfdbdat = griddata((a, b), dfdb, (adat, bdat), method=method)
    magdat = griddata((a, b), mag, (adat, bdat), method=method)

    # close the file if we opened it
    if opened_here:
        del pf

    # return in the proper order
    return adat, bdat, dfdadat, dfdbdat, magdat


ray_cache = {}

def lineout(p1, p2, field, pf, **kwargs):
    """Grabs a line out (ray) of a certain field (parameter, data) from the 
    flash output file pf.  

    This function requires both scipy and yt and returns x, y, z and value 
    data that is suitable for plotting with matplotlib.plot().
    
    Parameters
    ----------
    p1 : three-tuple of floats
        The first point in the line-out.
    p2 : three-tuple of floats
        The second point in the line-out.
    field : str
        A field to retrieve, e.g. 'dens'.
    pf : str or yt.data_objects.static_output
        The checkpoint or plot file, may be either a path to the file 
        or the yt object returned by pf.mods.load().
    kwargs : optional
        All other keyword arguments are passes to the ray() method
        on the pf.h object from yt.

    Returns
    -------
    x : 1D numpy array
        Interpolated x data.
    y : 1D numpy array
        Interpolated y data.
    z : 1D numpy array
        Interpolated z data.
    v : 1D numpy array
        Interpolated field values along ray.

    """
    # open the file if we got a path.
    opened_here = False
    if isinstance(pf, basestring):
        pf = yt.mods.load(pf)
        opened_here = True

    # make hashable
    p1 = tuple(p1)
    p2 = tuple(p2)
    kwargs = _to_hashable(kwargs)

    # get the ray object
    ray_key = (pf, p1, p2) + tuple(sorted(kwargs.items()))
    if ray_key in ray_cache:
        amr_ray = ray_cache[ray_key]
    else:
        amr_ray = pf.h.ray(p1, p2, field, **kwargs)
        ray_cache[ray_key] = amr_ray

    # get points
    x = amr_ray['x']
    y = amr_ray['y']
    z = amr_ray['z']
    v = amr_ray[field]

    # close the file if we opened it
    if opened_here:
        del pf

    # return in the proper order
    return x, y, z, v


def shock_on_lineout(p1, p2, field, pf, threshold=1e-6, min_threshold=1e-36, **kwargs):
    """Finds the shock of a certain field (parameter, data) along a line out (ray) 
    in a flash output file pf.  Currently this assumes that p1 is closer to the 
    center of the shock than p2.  Moreover the geometry from the FLASH simulation
    must be cartesian or 2D cylindrical (rz).

    Parameters
    ----------
    p1 : three-tuple of floats
        The first point in the line-out.
    p2 : three-tuple of floats
        The second point in the line-out.
    field : str
        A field to retrieve, e.g. 'dens'.
    pf : str or yt.data_objects.static_output
        The checkpoint or plot file, may be either a path to the file 
        or the yt object returned by pf.mods.load().
    threshold : float, optional
        The value above which gradients must exceed to be considered a
        shock.  Required for ignoring low-level noise.
    min_threshold : float, optional
        If a shock is not found at a given threshold level, the threshold
        value is reduced by an order of magnitude until a shock is found.
        This continues until a minimum threshold is reached.
    kwargs : optional
        All other keyword arguments are passes to the ray() method
        on the pf.h object from yt.

    Returns
    -------
    shock_p : three-tuple of floats
        Shock position.
    shock_v : field type (float)
        Shock peak value.

    See Also
    --------
    flash.output.lineout : 
        Used to get the lineout from two points and a data file.

    flash.analysis.shock_detect : 
        Used to find the shock itself.

    """
    x, y, z, v = lineout(p1, p2, field, pf)
    r = np.sqrt(x**2 + y**2 + z**2)
    sort_index = rad.argsort()
    r = r[sort_index]
    v = v[sort_index]
    shock_r, shock_v, shock_i = analysis.shock_detect(r, v, threshold=threshold, 
                                                      min_threshold=min_threshold)
    shock_p = (x[sort_index][shock_i], y[sort_index][shock_i], z[sort_index][shock_i])
    return shock_p, shock_v


def load_laser_dat(filename):
    """Reads in a LaserEnergyProfile.dat file as a numpy structured array.

    Parameters
    ----------
    filename : str or file-like
        The string path to the file or an object which implements the file
        interface

    Returns
    -------
    dat : ndarray
        This structured array has the following dtype and units:

        ========== ======= ========
        Field      dtype   units
        ========== ======= ========
        Step       int32   unitless
        Time       float64 s     
        dt         float64 s
        Energy In  float64 erg
        Energy Out float64 erg
        dE In      float64 erg
        dE Out     float64 erg
        ========== ======= ========

    """
    datdt = np.dtype([('step', np.int32),
                      ('t', np.float64),
                      ('dt', np.float64),
                      ('E_in', np.float64),
                      ('E_out', np.float64),
                      ('dE_in', np.float64),
                      ('dE_out', np.float64),
                      ])
    dat = np.loadtxt(filename, dtype=datdt, skiprows=2)
    return dat
