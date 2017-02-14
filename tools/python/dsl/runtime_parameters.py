"""Handles FLASH's runtime parameter file domain-specific language.  Since this language is 
basically a key-value store, we translate such files to/from Python dictionaries.  An AST 
representation is not needed."""

import re 
from itertools import izip_longest, groupby

# pre-compiled regexes
_par_line_pat = ("^\s*(\w+)\s*=\s*"
                 """([^'"][^#\s]*|'{1,3}?[^'\n]*?'{1,3}?|"{1,3}?[^"\n]*?"{1,3}?)"""
                 )
_par_line_reg = re.compile(_par_line_pat, re.M)

_indexed_key_reg = re.compile("(\w+)_(\d+)$")

_par_int_reg = re.compile("^[+-]?\d+$") 
_par_str_reg = re.compile("""^(".*"|'.*')$""") 
_par_bool_reg = re.compile("\.([Tt][Rr][Uu][Ee]|[Ff][Aa][Ll][Ss][Ee])\.") 
_par_float_reg = re.compile("^[+-]?(\d+\.?|\.?\d+)\d*[Ee]?[+-]?\d*$")

_py_match_convert = [
    (_par_str_reg, lambda v: v[1:-1]),
    (_par_int_reg, int),
    (_par_float_reg, float),
    (_par_bool_reg, lambda v: '.true.' == v.lower()),
    ]


def _grouper(n, iterable, fillvalue=None):
    """From itertools; grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"""
    args = [iter(iterable)] * n
    return list(izip_longest(fillvalue=fillvalue, *args))


def _flat_reshape(iterable, shape):
    vals = list(iterable)
    ndims = len(shape)

    nelems = 1
    for i in range(ndims):
        nelems *= shape[i]
    assert nelems == len(vals)

    for i in range(ndims - 1, 0, -1):
        vals = _grouper(shape[i], vals)

    return vals


def _nest_flat(iterable):
    """Takes a flat list whose elements are tuples
    of (index, value) items and returns a nested list
    of values, placed by the index."""
    # Ensure dense at this level
    i00 = set([i[0][0] for i in iterable])
    min_index = min(i00)
    max_index = max(i00)
    if len(i00) != max_index - min_index + 1:
        raise IndexError("iterable not dense at this level of nesting")

    # return final state
    if all([1 == len(i[0]) for i in iterable]):
        return [i[1] for i in iterable]

    # Nest by a single level
    nested = [sorted(b for b in a[1]) for a in groupby(iterable, lambda x: x[0][0])]

    # Remove an index level
    reduced_nested = []
    for nest in nested:
        reduced = [(n[0][1:], n[1]) for n in nest]
        reduced = _nest_flat(reduced)
        reduced_nested.append(reduced)

    return reduced_nested


def py_form(pardict):
    """Takes a dictionary of runtime parameters and returns a new dictionary
    whose values are in a canonically Python form."""
    pydict = {}
    indexed_pars = {}

    # roll through items, converting everything
    for key, value in pardict.items():
        # get base key and index, walks backwards from the end of the key str
        key_ind = (key, tuple())
        m = _indexed_key_reg.match(key)
        while m:
            basekey, ind = m.groups()
            key_ind = (basekey, (int(ind),) + key_ind[1])
            m = _indexed_key_reg.match(basekey)

        # convert value
        for regex, converter in _py_match_convert:
            if regex.match(value):
                pyvalue = converter(value)
                break
        else:
            raise ValueError("could not convert value in item: {0} = {1}".format(key, value))

        # Add to appropriate dictionary
        pykey, pyind = key_ind
        if 0 == len(pyind):
            pydict[pykey] = pyvalue
        else:
            if pykey not in indexed_pars:
                indexed_pars[pykey] = {}
            indexed_pars[pykey][pyind] = pyvalue

    # add indexed pairs as lists to pydict if possible
    for pykey in indexed_pars:
        pyitems = indexed_pars[pykey].items()
        pyitems.sort()

        # try to nest
        try:
            vals = _nest_flat(pyitems)
            isdense = True
        except IndexError:
            isdense = False

        # set the value(s)
        if isdense:
            # add to pydict as a nested list
            pydict[pykey] = vals
        else:
            # add to pydict as individual elements
            for ind, val in pyitems:
                itemkey = pykey + "_" + "_".join([str(i) for i in ind])
                pydict[itemkey] = val

    return pydict


def _isiterable(iterable):
    try:
        iter(iterable)
        b = True
    except TypeError:
        b = False
    return b


def _par_convert_iterable(key, value):
    flattened = {}
    for i, v in enumerate(value):
        flatkey = key + "_{0}".format(i + 1)
        for isfunc, extra_args, converter in _par_match_convert:
            if isfunc(v, *extra_args):
                flatvalues = converter(flatkey, v)
                break
        else:
            msg = "parameter type could not be inferred: {0} = {1}".format(flatkey, v)
            raise ValueError(msg)
        flattened.update(flatvalues)
    return flattened

_par_match_convert = [
    (isinstance, [basestring], lambda k, v: {k: '"{0}"'.format(v)}),
    (isinstance, [float], lambda k, v: {k: '{0:.8e}'.format(v)}),
    (isinstance, [bool], lambda k, v: {k: '.true.' if v else '.false.'}),  # Must be before int
    (isinstance, [int], lambda k, v: {k: str(v)}),
    (isinstance, [long], lambda k, v: {k: str(v)}),
    (_isiterable, [], _par_convert_iterable),
    ]

def par_form(pydict):
    """Takes a Python dictionary and returns a new dictionary
    whose values are in runtime parameters form."""
    pardict = {}
    for pykey, pyvalue in pydict.items():
        for isfunc, extra_args, converter in _par_match_convert:
            if isfunc(pyvalue, *extra_args):
                parvalues = converter(pykey, pyvalue)
                break
        else:
            msg = "parameter type could not be inferred: {0} = {1}".format(pykey, pyvalue)
            raise ValueError(msg)
        pardict.update(parvalues)
    return pardict


def load(parfile):
    """Loads a FLASH runtime parameter file, such as flash.par,
    and returns a Python dictionary."""
    # get the file handler, if not provided
    opened_here = False
    if isinstance(parfile, basestring):
        parfile = open(parfile)
        opened_here = True

    # read the file in.
    rawfile = parfile.read()
    if opened_here:
        parfile.close()

    # parse file into dictionary
    rawdict = dict(_par_line_reg.findall(rawfile))

    # transform to python-form
    pydict = py_form(rawdict)

    return pydict


def dump(pydict, parfile="flash.par", mode='w'):
    """Writes a FLASH runtime parameter file, such as flash.par,
    from a Python dictionary.  Such dictionaries mus be fairly 
    simple as the runtime parameters file specification is not 
    very sophisticatd."""
    # convert dict to par-form
    pardict = par_form(pydict)

    # convert to string
    paritems = pardict.items()
    paritems.sort()
    paritems = ["{0} = {1}".format(k, v) for k, v in paritems]
    rawfile = "\n".join(paritems)

    # get the file handler, if not provided
    opened_here = False
    if isinstance(parfile, basestring):
        parfile = open(parfile, mode)
        opened_here = True

    # write the file
    parfile.write(rawfile)
    if opened_here:
        parfile.close()

