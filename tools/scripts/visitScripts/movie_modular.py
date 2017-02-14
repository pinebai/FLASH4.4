# #!/usr/bin/env python
#!/usr/bin/python

#================================================================================
#  movie_modular.py
#================================================================================
#  This script makes a series of images from a time series of flash files on davinci.
#  The images can then be combined into a movie.
#
#  FUNCTIONS:
#    range_intersection(): finds intersection of ranges of values of a variable over all
#                          files in series
#    range_union(): finds union of ranges of values of a variable over all
#                   files in series
#    snap(): loops over files and makes images
#    set_plots_ops_1(): adds pseudocolor and vector plots, and operators, and sets their attributes
#    set_plots_ops_2(): adds contour plot, and operators, and sets their attributes
#    set_ref_level(): chooses the refinement level to plot
#================================================================================

import glob
import sys
import os
import subprocess
import os.path
# import optparse

#================================================================================
#================================================================================
#   TO TRY TO KEEP HEDGEHOG SIZE RIGHT FOR WHOLE SERIES
# VARY "N vectors" of vec attributes, START TO FINISH, FROM 23 TO 11
# VARY "Scale" of vec attributes, START TO FINISH, FROM 0.00015 TO 0.0003
# VARY "Zoom" of  from ~300 to ~250
# VARY "Focus" to between (0, 4.7e6, 170) in beginning, and (0, 4.5e7, -30000)
#================================================================================
#================================================================================


##################################################################################
## I HAVEN'T USED THESE FOR A WHILE....
##################################################################################
##................................................................................
## SESSION FILE:
#sessionfilename = "movietest_rpv1contour_1.session"
##................................................................................
## IMAGE FILE:
## imagefilename = "/Users/hudson/Pictures/VisIt/BubbleScript/OnLiberty/testimage"
#if_path = "/Users/hudson/Pictures/VisIt/BubbleScript/OnLiberty/"
#ifn_prefix = "rpv1slice_velhhogs_"
#
## PLOT & IMAGE FILES:
## step_start = 0
## step_end = 300
## step_incr = 10
##................................................................................
##################################################################################


#================================================================================
# HARD-CODED PARAMETERS THAT STAND IN FOR OPTPARSE ONES UNTIL OPTPARSE IS AVAILABLE
#
#................................................................................
# dbase_loc:   1) indicates where database is (and, therefore, where compute engine will be launched)
#    and 2) determines which command is used to read database names and some metadata
#    before calling visit methods
# client_loc:  1) indicates where client will run and 2) determines which command is
#    used to launch compute engine
# movie_id:    indicates which combo of plots & ops should be displayed
#................................................................................
dbase_loc = "davinci" # Location of database (and compute engine) # "local"
client_loc = "davinci" # Location of client # "local" "davinci"
#engine_loc = "davinci" # "davinciremote" ==> client local, compute engine on davinci
                            # "davincilocal" ==> client & compute engine on davinci
                            # "local" ==> client & compute engine local
                            # "teragrid", "ellipse"
movie_id = "hhogs_isosurf_clipped" # "sline_slice":  streamlines (and rpv1 contour) in an X=0 slice
                         # "hhogs_isosurf_clipped":  contour surface of rpv1
#                                                  + edge of that contour
#                                                  + thresholded hedgehog plot of velocity
#................................................................................
num_dims = "3d" # "2d", "3d"
commonrange = "True" # RANGE_INTERSECTION() IS CURRENTLY NOT OPTIONAL, SINCE IT ALSO GLEANS
                     # FROM FILE NAMES THE TIME STEP NUMBERS....

#................................................................................
# PLOT FILE ("pf" means "plotfile"):
# pf_path = "/project/projectdirs/incite13/analysis_cal/500m_data_extract/"
pf_path = "/project/projectdirs/incite13/hudson/Data/FromHPSS/"
# pf_path = "/project/projectdirs/incite13/hudson/Data/TEMP/"
# pf_path = "/Users/hudson/Data/3d/"
#................................................................................

rmt_machine_pfx = "davinci.nersc.gov:"
# pf_remotepath = "davinci.nersc.gov:/project/projectdirs/incite13/analysis_cal/500m_data_extract/"
# davinci_pfpath = "davinci.nersc.gov:/project/projectdirs/incite13/analysis_cal/500m_data_extract/"

#................................................................................
pfn_prefix = "500m_16o40_hdf5_plt_cnt_"
# pfn_prefix = "c_hdf5_plt_cnt_"
#................................................................................

#................................................................................
vname = "rpv1"
# vname = "dens"
#................................................................................

filelist = []

# OTHER
# h5path = "/usr/common/usg/hdf5/1.6.4/serial/bin/"
h5path = "/usr/common/homes/b/brugger/visit/hdf5/1.6.5/altix_gcc_3.3.3/bin/"
# h5path = "/Users/hudson/Software/hdf5_fromsource/hdf5-1.6.5/hdf5/bin/"

viewscope_flag = 1
camvec = (1, 0.5, -0.75) # (1, 0, 0) # (1, 0, 0.5)
#--------------------------------------------------------------------------------
zoom = 1
zoom_start = 5 # 647 # 25 # 1000 # 300 # 200
zoom_stop = 5 # 133 # 25 # 250 # 250 # 350
zoom_incr = 0
#--------------------------------------------------------------------------------
focus_start = (0, 16e6, 0) # (0, 5e7, 170) # (0, 4.7e6, 170) # (0, 2.0e7, 0.0)
focus_stop = (0, 16e6, 0) # (0, 2e7, 170) # (0, 4.5e7, -30000) # (0, 3.5e7, 0.0)
focus_incr = (0.0, 0.0, 0.0)
focus = []
#--------------------------------------------------------------------------------
upvec = []
upvec_plus_y = [0.0, 1.0, 0.0]
upvec_minus_z = [0.0, 0.0, -1.0]
#--------------------------------------------------------------------------------
state_incr = 1 # The increment to the next state to save as image in the slider loop
#--------------------------------------------------------------------------------
pan = (0.0, 0.0)
pan_start = (0.0, 0.0) # (1.71e-5, -0.00087)
pan_stop = (0.0, 0.0) # (-1.27e-5, -0.0038)
pan_incr = (0.0, 0.0)
#--------------------------------------------------------------------------------
numvecs = 0
numvecs_start = 23 # 23
numvecs_stop = 9 # 11
numvecs_incr = 0
vecscale = 1.0
vecscale_start = 0.0002 # 0.00015
vecscale_stop = 0.0004 # 0.0003
vecscale_incr = 0.0
#--------------------------------------------------------------------------------
# FOR STREAMLINE
#--------------------------------------------------------------------------------
radius_factor = 0.8
sline_StepLength = 0.01
sline_MaxSteps = 200.0
sline_PointDensity = 15
#sline_ColorBy = speed # Default
sline_ShowStart = 0
sline_BoxX = [-2e6, 1.5e7]
#sline_DisplayAs = Tubes

clip_Plane1Origin = (1e5, 0.0, 0.0) # For the slab of streamlines
clip_Plane1Normal = (1, 0, 0)
clip_Plane2Origin = (1e5, 0.0, 0.0)
clip_Plane2Normal = (0, 1, 0)
clip1_Plane1Origin = (0.0, 0.0, 0.0) # For the two Contour plots in set_plots_8()
clip1_Plane1Normal = (1, 0, 0)
clip1_Plane2Origin = (0.0, 0.0, 0.0)
clip1_Plane2Normal = (0, 1, 0)
clip2_Plane1Origin = (0.0, 0.0, 0.0) # For the Vector plot in set_plots_8()
clip2_Plane1Normal = (1, 0, -1) # (1, 0, 0)
clip2_Plane2Origin = (0.0, 0.0, 0.0)
clip2_Plane2Normal = (0, 1, 0)

resample_StartEndX = [-2e6, 1.5e7]
resample_SamplesXYZ = [40, 20, 40]
#--------------------------------------------------------------------------------
#................................................................................
#................................................................................
# For cubic spline interpolation of vector stride and hedeghog scale between
#   keyframe values of those parameters
#   (head and width of hhog abide)
#   values need to be floats for EvalCubicSpline()
#
# EvalCubicSpline(t, weights, values) -> f(t)
#   t:        A floating point value in the range [0., 1.] that represents the distance from
#             the first control point to the last control point.
#   weights:  A tuple of N floating point values in the range [0., 1.] that represent how
#             far along in parameterized space, the values will be located.
#   values:   A tuple of N objects to be blended. Any objects that can be used in
#             arithmetic expressions can be passed.
#
# Prefix "kf_" means "keyframe"
# "timesteps" contains, as floats, the step numbers of ALL files processed
# For interpolation of vector attributes:
#    "kf_timesteps", "kf_vecstride" and "kf_hhogscale" must have same # elements
#    "kf_vecstride" and "kf_hhogscale" will act as "values" tuples in 2 calls to EvalCubicSpline()
# For interpolation of streamline attributes:
#    The following must have same # elements.  All but "kf_timesteps" will act as "values" tuples in 2 calls to EvalCubicSpline().
#    "kf_timesteps"
#    "kf_SlineRadius"
#    "kf_ResampleStartY"
#    "kf_ResampleEndY"
#    "kf_ResampleStartZ"
#    "kf_ResampleEndZ"
#    "kf_SlineBoxYLo"
#    "kf_SlineBoxYHi"
#    "kf_SlineBoxZLo"
#    "kf_SlineBoxZHi"
# "timesteps" contains, as floats, the step numbers of ALL files processed
kf_timesteps = (1.0, 49.0, 101.0, 149.0, 201.0, 249.0, 305.0)
#................................................................................
kf_vecstride = (5.0, 14.0, 25.0, 40.0, 75.0, 215.0, 305.0) # Values tediously, manually determined
kf_hhogscale = (0.0001, 0.0002, 0.0003, 0.00036, 0.00065, 0.001, 0.0013) # Values tediously, manually determined
#................................................................................
kf_SlineRadius = (3.5e4, 5e4, 1e5, 1.25e5, 1.75e5, 2e5, 2.25e5)
kf_ResampleStartY = (2.5e6, 2.5e6, 3e6, 7e6, 1.5e7, 2.5e7, 3.8e7)
kf_ResampleEndY = (5.5e6, 7.5e6, 1.1e7, 1.6e7, 2.5e7, 4e7, 6e7)
kf_ResampleStartZ = (-1.5e6, -2.5e6, -4e6, -5e6, -7e6, -8e6, -1.5e7)
kf_ResampleEndZ = (1.5e6, 2.5e6, 4e6, 5e6, 7e6, 8e6, 1.5e7)
kf_SlineBoxYLo = (2.5e6, 2.5e6, 3e6, 7e6, 1.5e7, 2.5e7, 3.8e7)
kf_SlineBoxYHi = (5.5e6, 7.5e6, 1.1e7, 1.6e7, 2.5e7, 4e7, 6e7)
kf_SlineBoxZLo = (-1.5e6, -2.5e6, -4e6, -5e6, -7e6, -8e6, -1.5e7)
kf_SlineBoxZHi = (1.5e6, 2.5e6, 4e6, 5e6, 7e6, 8e6, 1.5e7)
#................................................................................

num_keyframes = len(kf_timesteps)
num_movieframes = 77
first_movieframe = 1.0
last_movieframe = 305.0
timesteps = []
fnames = [] # Global only for diagnosing allocation crash after beaucoup successfully processing time steps
#................................................................................
#................................................................................
#--------------------------------------------------------------------------------
#streamplaneorigin = (0, 0, 0)
#streamplaneorigin_start = (0, 4.5e6, 0)
#streamplaneorigin_stop = (0, 4.5e7, 0)
#streamplaneorigin_incr = (0, 0, 0)
streamplaneorigin = [0.0, 0.0, 0.0]
streamplaneorigin_start = [0.0, 4.5e6, 0.0] # [0.0, 3.0e7, 0.0]
streamplaneorigin_stop = [0.0, 4.5e7, 0.0] # [0.0, 3.0e7, 0.0]
streamplaneorigin_incr = [0.0, 0.0, 0.0]
streamplaneradius = 1.0
streamplaneradius_start = 2.0e6 # 2.0e6 # 1.0e7
streamplaneradius_stop = 1.0e7 # 1.5e7 # 1.0e7
streamplaneradius_incr = 0.0
#--------------------------------------------------------------------------------
#================================================================================

minima = []
maxima = []
lowmax = highmin = 0.0
lowmin = highmax = 0.0

#================================================================================
# Use visit Query() to determine files' mins and maxes
def queryminmax():
    print "About to try looping over database for mins and maxes"
#    for state in range(TimeSliderGetNStates()):
#        SetTimeSliderState(state)
#        DrawPlots()
#        print "State ", state
#        Query("MinMax")
#        print "Just queried"
#        print GetQueryOutputString()
#        print "Just printed query string"
#        print GetQueryOutputValue()

#================================================================================
# Determine the value range common to all files in a series by reading files w/ h5dump
def range_intersection():
    global minima, maxima
    global lowmax, highmin

    # ---------------------------------------------------------------------------
    #  BUILD LIST OF NAMES OF FLASH FILES
#    pf_remotepath = pf_remoteprefix + pf_path

#    if engine_loc == "davinci":
    if dbase_loc == "davinci":
        if client_loc == "davinci":
            lscmd = "ls " + pf_path + pfn_prefix + "*"
            print "[range_intersection()] - cmd to list LOCAL files", lscmd
        elif client_loc == "local":
            lscmd = "ssh davinci.nersc.gov -l rhudson \"ls " + pf_path + pfn_prefix + "*\""
            print "[range_intersection()] - cmd to list REMOTE files", lscmd
#    elif engine_loc == "local":
    elif dbase_loc == "local":
        lscmd = "ls " + pf_path + pfn_prefix + "*"
        print "[range_intersection()] - cmd to list LOCAL files", lscmd

    p1 = subprocess.Popen(lscmd, shell=True, stdout=subprocess.PIPE)
    filelist = p1.stdout.read()
    sl_filelist = filelist.splitlines()
    print "----------------------------------------------------------------------"
    print "[range_intersection()] - file names with path: ", sl_filelist
    print "----------------------------------------------------------------------"

    for thispath in sl_filelist:
        fname = os.path.basename(thispath)
        fnames.append(fname)
    print "----------------------------------------------------------------------"
    print "[range_intersection()] - file names without path: ", fnames
    print "----------------------------------------------------------------------"

    for fname in fnames:
        fname = fname.split('_')[-1]
        timesteps.append(float(fname))
    print "----------------------------------------------------------------------"
    print "[range_intersection()] - time steps: ", timesteps
    print "[range_intersection()] - number of time steps: ", len(timesteps)
    print "----------------------------------------------------------------------"
#    sys.exit()

    numfiles = len(sl_filelist)
    print "[range_intersection()] - Files in series: ", filelist
    print "[range_intersection()] - Num files in series: ", numfiles

    # ---------------------------------------------------------------------------
    #  IF MIN & MAX ARE NOT STORED IN FILES, EXIT (test existence of "mininum" attribute in first file)
    min_attributepath = "/" + vname + "/" + "minimum"
#    min_h5cmd = "/usr/common/usg/hdf5/1.6.4/serial/bin/h5dump -a " + "%s %s"%(min_attributepath, sl_filelist[0])
    min_h5cmd = h5path + "/h5dump -a " + "%s %s"%(min_attributepath, sl_filelist[0])

#    if engine_loc == "davinci":
##        rmt_min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
#        min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
##        ret = subprocess.call(rmt_min_h5cmd, shell=True)
##    elif engine_loc == "local":
##        ret = subprocess.call(min_h5cmd, shell=True)
    if dbase_loc == "davinci" and client_loc == "local":
        min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
    # otherwise, defaults to local path

    ret = subprocess.call(min_h5cmd, shell=True)
    if ret != 0:
        print "File %s contains no attribute \"minimum\" (or, I conclude, \"maximum\")", sl_filelist[0]
        print "Exiting..."
        sys.exit()

    # ---------------------------------------------------------------------------
    #  GET MINIMA AND MAXIMA FROM FLASH FILES
    max_attributepath = "/" + vname + "/" + "maximum"
    for l in range(1, len(sl_filelist)):  #  Skip first file, whose max and min may == 0.0
#        min_h5cmd = "/usr/common/usg/hdf5/1.6.4/serial/bin/h5dump -a " + "%s %s"%(min_attributepath, sl_filelist[l])
#        max_h5cmd = "/usr/common/usg/hdf5/1.6.4/serial/bin/h5dump -a " + "%s %s"%(max_attributepath, sl_filelist[l])
        min_h5cmd = h5path + "/h5dump -a " + "%s %s"%(min_attributepath, sl_filelist[l])
        max_h5cmd = h5path + "/h5dump -a " + "%s %s"%(max_attributepath, sl_filelist[l])
#        rmt_min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
#        rmt_max_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + max_h5cmd + "\""
#        if engine_loc == "davinci":
        if dbase_loc == "davinci" and client_loc == "local":
            min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
            max_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + max_h5cmd + "\""

        #  MINIMUM
#        pmin = subprocess.Popen(rmt_min_h5cmd, shell=True, stdout=subprocess.PIPE)
        pmin = subprocess.Popen(min_h5cmd, shell=True, stdout=subprocess.PIPE)
        minout = pmin.stdout.read()
        sl_minout = minout.splitlines()
        for i in range(len(sl_minout)):
            s_sl_minout = sl_minout[i].split()
            if s_sl_minout[0] == "(0):":
                fmin = float(s_sl_minout[1])
                minima.append(fmin)

        #  MAXIMUM
#        pmax = subprocess.Popen(rmt_max_h5cmd, shell=True, stdout=subprocess.PIPE)
        pmax = subprocess.Popen(max_h5cmd, shell=True, stdout=subprocess.PIPE)
        maxout = pmax.stdout.read()
        sl_maxout = maxout.splitlines()
        for i in range(len(sl_maxout)):
            s_sl_maxout = sl_maxout[i].split()
            if s_sl_maxout[0] == "(0):":
                fmax = float(s_sl_maxout[1])
                maxima.append(fmax)

    print "[range_intersection()] - All minima: ", minima
    print "[range_intersection()] - All maxima: ", maxima

    lowmax = min(maxima)
    highmin = max(minima)
    print "[range_intersection()] - Lowest maximum: ", lowmax
    print "[range_intersection()] - Highest minimum: ", highmin

#================================================================================
# Determine...whatever
# THIS IS NOT DONE, YET; NEEDS TO BE UPGRADED IN THE SAME MANNER AS RANGE_INTERSECTION
def range_union():
    global minima, maxima
    global lowmax, highmin

    #  GET NAMES OF REMOTE FLASH FILES
#    pf_remotepath = pf_remoteprefix + pf_path
    lscmd = "ssh davinci.nersc.gov -l rhudson \"ls " + pf_path + pfn_prefix + "*\""
    p1 = subprocess.Popen(lscmd, shell=True, stdout=subprocess.PIPE)
    filelist = p1.stdout.read()
    sl_filelist = filelist.splitlines()
    numfiles = len(sl_filelist)
    print "Num files in series: ", numfiles

    #  IF MIN & MAX ARE NOT STORED IN FILES, EXIT (test existence of "mininum" attribute in first file)
    min_attributepath = "/" + vname + "/" + "minimum"
    min_h5cmd = "/usr/common/usg/hdf5/1.6.4/serial/bin/h5dump -a " + "%s %s"%(min_attributepath, sl_filelist[0])
    rmt_min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
    ret = subprocess.call(rmt_min_h5cmd, shell=True)
    if ret != 0:
        print "File %s contains no attribute \"minimum\" (or, I conclude, \"maximum\")", sl_filelist[0]
        print "Exiting..."
        sys.exit()

    #  GET MINIMA AND MAXIMA FROM REMOTE FLASH FILES
    max_attributepath = "/" + vname + "/" + "maximum"
    for l in range(1, len(sl_filelist)):  #  Skip first file, whose max and min may == 0.0
        min_h5cmd = "/usr/common/usg/hdf5/1.6.4/serial/bin/h5dump -a " + "%s %s"%(min_attributepath, sl_filelist[l])
        max_h5cmd = "/usr/common/usg/hdf5/1.6.4/serial/bin/h5dump -a " + "%s %s"%(max_attributepath, sl_filelist[l])
        rmt_min_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + min_h5cmd + "\""
        rmt_max_h5cmd = "ssh davinci.nersc.gov -l rhudson \"" + max_h5cmd + "\""

        #  MINIMUM
        pmin = subprocess.Popen(rmt_min_h5cmd, shell=True, stdout=subprocess.PIPE)
        minout = pmin.stdout.read()
        sl_minout = minout.splitlines()
        for i in range(len(sl_minout)):
            s_sl_minout = sl_minout[i].split()
            if s_sl_minout[0] == "(0):":
                fmin = float(s_sl_minout[1])
                minima.append(fmin)

        #  MAXIMUM
        pmax = subprocess.Popen(rmt_max_h5cmd, shell=True, stdout=subprocess.PIPE)
        maxout = pmax.stdout.read()
        sl_maxout = maxout.splitlines()
        for i in range(len(sl_maxout)):
            s_sl_maxout = sl_maxout[i].split()
            if s_sl_maxout[0] == "(0):":
                fmax = float(s_sl_maxout[1])
                maxima.append(fmax)

    print "[range_union()] - All minima: ", minima
    print "[range_union()] - All maxima: ", maxima

    lowmin = min(minima)
    highmax = max(maxima)
    print "[range_union()] - Lowest minimum: ", lowmin
    print "[range_union()] - Highest maximum: ", highmax


#================================================================================
# Establish view bounds as either entire domain or merely the plots
def set_viewscope():
    # ---------------------------------------------------------------------------
    # Select scope of view (bounds 0) entire domain or 1) only plotted data)
    if viewscope_flag == 0:
        SetViewExtentsType("original") # Show entire data dmoain
    elif viewscope_flag == 1:
        SetViewExtentsType("actual") # Crops view tightly about the plot


#================================================================================
# Set up plots and operators for pseudocolor & vectors
def set_plots_ops_1():
#    global vec_atts
    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Pseudocolor", visit_vname)
#    AddPlot("Vector", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
#    l_atts.type = 1
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)

    # Pseudocolor ---------------------------
    pc_atts = PseudocolorAttributes()
    pc_atts.centering = pc_atts.Nodal # center the values at cell corners
    pc_atts.opacity = 1.0
#    pc_atts.limitsMode = 1 # limits a function of 0) original data or 1) current plot
    if commonrange == "True":
        pc_atts.min = highmin # lower bound of range that's colormapped
        pc_atts.max = lowmax # upper bound of range that's colormapped
        pc_atts.minFlag = 1 # turn on
        pc_atts.maxFlag = 1 # turn on
    SetPlotOptions(pc_atts)

    # Vector ---------------------------------
#    vec_atts = VectorAttributes()
#    vec_atts.useStride = 1
#    vec_atts.stride = 12
#    vec_atts.scale = 0.00025
#    vec_atts.scaleByMagnitude = 0
#    vec_atts.headSize = 0.22
#    vec_atts.lineStem = 1
#    vec_atts.highQuality = 1
#    vec_atts.stemWidth = 0.15
#    SetPlotOptions(vec_atts)

    # Set up operator(s)
    if num_dims == "3d":
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
#        s_atts.originPoint = (0, 0, 1.5e7)
#        s_atts.originType=s_atts.Point
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
#        xform_atts.translateZ = 2e7
#        xform_atts.doTranslate = 1
        SetOperatorOptions(xform_atts)

#    SetActivePlots(0) # Make only Pseudocolor plot active for Isosurface and Tube operators
#    AddOperator("Isosurface")
#    iso_atts = IsosurfaceAttributes()
#    iso_atts.contourMethod = 1
#    iso_atts.contourValue = 0.5
#    SetOperatorOptions(iso_atts)
#    AddOperator("Tube")
#    tube_atts = TubeAttributes()
#    tube_atts.width = 120000 # 75000
#    tube_atts.fineness = 8
#    tube_atts.capping = 1
#    SetOperatorOptions(tube_atts)
#    SetActivePlots((0, 1)) # Make both plots active

##    SetActivePlots(1) # Make only Vector plot active for Threshold operator
##   THRESHOLD OPERATOR
#    AddOperator("Threshold")
#    t_atts = ThresholdAttributes()
##    t_atts.variable = "mesh_blockandlevel/rpv1"
#    t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
#    print "----------------------------------------------------------------------"
#    print "[set_plots_ops_1()] - t_atts.defaultVarName", t_atts.defaultVarName
#    print "----------------------------------------------------------------------"
##    t_atts.lbound = 0.45
#    t_atts.lowerBounds = 0.45
#    print "----------------------------------------------------------------------"
#    print "[set_plots_ops_1()] - t_atts.lowerBounds", t_atts.lowerBounds
#    print "----------------------------------------------------------------------"
#    SetOperatorOptions(t_atts)

#    SetActivePlots((0, 1)) # Make both plots active

#    SetActivePlots(0) # Make only Pseudocolor plot active for removing Threshold operator
#    RemoveOperator(4)
#    SetActivePlots((0, 1)) # Make both plots active

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window


#================================================================================
# Set up plots and operators for contour
def set_plots_ops_2():
    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Contour", visit_vname)
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)

    # Contour --------------------------------
    c_atts = ContourAttributes()
    c_atts.contourMethod = c_atts.Value # center the values at cell corners
    c_atts.contourValue = 0.5
    c_atts.lineWidth = 3
#    c_atts.limitsMode = 1 # limits a function of 0) original data or 1) current plot

#    if commonrange == "True":
#        c_atts.min = highmin # lower bound of range that's colormapped
#        c_atts.max = lowmax # upper bound of range that's colormapped
#        c_atts.minFlag = 1 # turn on
#        c_atts.maxFlag = 1 # turn on
    SetPlotOptions(c_atts)

    # Set up operator(s)
    if num_dims == "3d":
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
        SetOperatorOptions(xform_atts)

#   THRESHOLD OPERATOR
#    AddOperator("Threshold")
#    t_atts = ThresholdAttributes()
#    t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
#    t_atts.lowerBounds = 0.45
#    SetOperatorOptions(t_atts)

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window


#================================================================================
# Set up vector plot and operators
def set_plots_ops_3():
    global numvecs, numvecs_start, vecscale, vecscale_start
#    global vec_atts

    # ---------------------------------------------------------------------------
    # Set up plot(s)
#    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Vector", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
#    l_atts.type = 1
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)

    # Vector ---------------------------------
#    vec_atts = VectorAttributes()
    # -----------------
#    vec_atts.useStride = 1
#    vec_atts.stride = 12
#    vec_atts.scale = 0.00025
#    vec_atts.scaleByMagnitude = 0
#    vec_atts.headSize = 0.22
#    vec_atts.lineStem = 1
#    vec_atts.highQuality = 1
#    vec_atts.stemWidth = 0.15
    # -----------------
    vec_atts.useStride = 0
    numvecs = numvecs_start
    vec_atts.nVectors = numvecs
    vecscale = vecscale_start
    vec_atts.scale = vecscale
    vec_atts.scaleByMagnitude = 0
    vec_atts.headSize = 0.22
    vec_atts.lineStem = 1
    vec_atts.highQuality = 1
    vec_atts.stemWidth = 0.15
    # -----------------
    SetPlotOptions(vec_atts)

    # Set up operator(s)
    if num_dims == "3d":
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
#        s_atts.originPoint = (0, 0, 1.5e7)
#        s_atts.originType=s_atts.Point
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
#        xform_atts.translateZ = 2e7
#        xform_atts.doTranslate = 1
        SetOperatorOptions(xform_atts)

#   THRESHOLD OPERATOR
    AddOperator("Threshold")
    t_atts = ThresholdAttributes()
#    t_atts.variable = "mesh_blockandlevel/rpv1"
    t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
    print "----------------------------------------------------------------------"
    print "[set_plots_ops_3()] - t_atts.defaultVarName", t_atts.defaultVarName
    print "----------------------------------------------------------------------"
#    t_atts.lbound = 0.45
    t_atts.lowerBounds = 0.45
    print "----------------------------------------------------------------------"
    print "[set_plots_ops_3()] - t_atts.lowerBounds", t_atts.lowerBounds
    print "----------------------------------------------------------------------"
    SetOperatorOptions(t_atts)

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window


#================================================================================
# Set up plots and operators for streamlines
def set_plots_ops_4():
    # ---------------------------------------------------------------------------
    # Set up plot(s)
#    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Streamline", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)

    # Streamline --------------------------------
    sl_atts = StreamlineAttributes()
    sl_atts.sourceType = sl_atts.SpecifiedPlane # center the values at cell corners
    sl_atts.stepLength = 0.1
    sl_atts.maxTime = 25.0
    sl_atts.planeOrigin = (0.0, 2.0e7, 0.0)
    sl_atts.planeNormal = (1.0, 0.0, 0.0)
    sl_atts.planeUpAxis = (0.0, 1.0, 0.0)
    sl_atts.planeRadius = 8.0e6
    sl_atts.pointDensity = 15 # 10
    sl_atts.displayMethod = sl_atts.Tubes
    sl_atts.showStart = 1
    sl_atts.radius = 1.0e5
    sl_atts.coloringMethod = sl_atts.ColorBySpeed
    SetPlotOptions(sl_atts)

    # Set up operator(s)
#    if num_dims == "3d":
#        AddOperator("Slice")
#        s_atts = SliceAttributes()
#        s_atts.project2d = 0
#        s_atts.normal = (1,0,0)
#        s_atts.upAxis = (0,1,0)
#        SetOperatorOptions(s_atts)
#        AddOperator("Transform")
#        xform_atts = TransformAttributes()
#        xform_atts.rotateAxis = (1, 0, 0)
#        xform_atts.rotateAmount = -90
#        xform_atts.doRotate = 1
#        SetOperatorOptions(xform_atts)

#   THRESHOLD OPERATOR
#    AddOperator("Threshold")
#    t_atts = ThresholdAttributes()
#    t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
#    t_atts.lowerBounds = 0.45
#    SetOperatorOptions(t_atts)

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window


#================================================================================
# Set up plots and operators for streamlines and contour slice
def set_plots_ops_5():
#    global sl_atts
    global streamplaneorigin, streamplaneorigin_start, streamplaneradius, streamplaneradius_start

    print "[set_plots_ops_5()] - 0"

    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
#    AddPlot("Contour", visit_vname)
    AddPlot("Contour", "mag_velocity_3d")
    AddPlot("Streamline", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations
#    print "[set_plots_ops_5()] - 1"

    l_atts = LightAttributes()
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)
#    print "[set_plots_ops_5()] - 2"

    # Contour --------------------------------
    c_atts = ContourAttributes()
#    c_atts.contourMethod = c_atts.Value # center the values at cell corners
#    c_atts.contourValue = 0.8 # 0.5
    c_atts.contourMethod = c_atts.Percent # center the values at cell corners
    c_atts.contourPercent = 75.0
    c_atts.lineWidth = 3
    SetPlotOptions(c_atts)

    # Streamline --------------------------------
    sl_atts = StreamlineAttributes()
    sl_atts.sourceType = sl_atts.SpecifiedPlane # center the values at cell corners
    sl_atts.stepLength = 0.05 # 0.01
    sl_atts.maxTime = 25.0
    #
    streamplaneorigin = streamplaneorigin_start
    print "-----------------------------------------------------------------------"
    print "[set_plots_ops_5()] - streamplaneorigin: ", streamplaneorigin
    print "-----------------------------------------------------------------------"
    sl_atts.planeOrigin = (streamplaneorigin[0], streamplaneorigin[1], streamplaneorigin[2])
    #
    sl_atts.planeNormal = (1.0, 0.0, 0.0)
    sl_atts.planeUpAxis = (0.0, 1.0, 0.0)
    #
    streamplaneradius = streamplaneradius_start
    print "-----------------------------------------------------------------------"
    print "[set_plots_ops_5()] - streamplaneradius: ", streamplaneradius
    print "-----------------------------------------------------------------------"
    sl_atts.planeRadius = streamplaneradius
    #
    sl_atts.pointDensity = 4 # 10 # 15
    sl_atts.displayMethod = sl_atts.Tubes
    sl_atts.showStart = 1
    sl_atts.radius = 1.0e5
    sl_atts.coloringMethod = sl_atts.ColorBySpeed
    SetPlotOptions(sl_atts)
#    print "[set_plots_ops_5()] - 3"

    # Set up operator(s)
    if num_dims == "3d":
        SetActivePlots(0) # Make only Contour plot active for Slice op
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)
        SetActivePlots((0, 1)) # Make both plots active
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
        SetOperatorOptions(xform_atts)

#    print "[set_plots_ops_5()] - 4"

##   THRESHOLD OPERATOR
#    SetActivePlots(1) # Make only Streamline plot active for Threshold op
#    AddOperator("Threshold")
#    t_atts = ThresholdAttributes()
#    t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
#    t_atts.lowerBounds = 0.45
##    t_atts.upperBounds = 0.55
#    SetOperatorOptions(t_atts)
#    SetActivePlots((0, 1)) # Make both plots active

    DrawPlots()
#    print "[set_plots_ops_5()] - 5"
    RedrawWindow() # (re)enable redrawing of vis window
#    print "[set_plots_ops_5()] - 6"

    return sl_atts


#================================================================================
# Set up plots and operators for contour & hedgehogs
def set_plots_ops_6():
    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Contour", visit_vname)
    AddPlot("Vector", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
    l_atts.direction = (1, 0, -1) # (-1, 0, -1)
    SetLight(0, l_atts)

    # Contour --------------------------------
    c_atts = ContourAttributes()
    c_atts.contourMethod = c_atts.Value # center the values at cell corners
    c_atts.contourValue = 0.5
    c_atts.lineWidth = 3
#    c_atts.limitsMode = 1 # limits a function of 0) original data or 1) current plot
#    if commonrange == "True":
#        c_atts.min = highmin # lower bound of range that's colormapped
#        c_atts.max = lowmax # upper bound of range that's colormapped
#        c_atts.minFlag = 1 # turn on
#        c_atts.maxFlag = 1 # turn on
    SetPlotOptions(c_atts)

    # Vector ---------------------------------
    vec_atts = VectorAttributes()
#    vec_atts.useStride = 0
#    numvecs = numvecs_start
#    vec_atts.nVectors = 70
    vec_atts.useStride = 1
    vec_atts.stride = int(kf_vecstride[0])
    vec_atts.scale = kf_hhogscale[0]
    vec_atts.scaleByMagnitude = 0
    vec_atts.headSize = 0.22
    vec_atts.lineStem = 0
    vec_atts.highQuality = 1
    vec_atts.stemWidth = 0.19
    # -----------------
    SetPlotOptions(vec_atts)

    # Set up operator(s)
    SetActivePlots((0, 1)) # Make both plots active (This seems needed here to make both plots active)
    if num_dims == "3d":
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
        SetOperatorOptions(xform_atts)

#   THRESHOLD OPERATOR
    SetActivePlots(1) # Make only Vector plot active for Threshold op
    AddOperator("Threshold")
    t_atts = ThresholdAttributes()
    t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
    t_atts.lowerBounds = 0.45
#    t_atts.upperBounds = 0.9999
    SetOperatorOptions(t_atts)
    SetActivePlots((0, 1)) # Make both plots active

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window

    return vec_atts


#================================================================================
# Set up plots and operators for contour & streamlines
def set_plots_ops_7():
    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Contour", visit_vname)
    AddPlot("Streamline", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)

    # Contour --------------------------------
    c_atts = ContourAttributes()
    c_atts.contourMethod = c_atts.Value # center the values at cell corners
    c_atts.contourValue = 0.5
    c_atts.lineWidth = 3
#    c_atts.limitsMode = 1 # limits a function of 0) original data or 1) current plot
#    if commonrange == "True":
#        c_atts.min = highmin # lower bound of range that's colormapped
#        c_atts.max = lowmax # upper bound of range that's colormapped
#        c_atts.minFlag = 1 # turn on
#        c_atts.maxFlag = 1 # turn on
    SetPlotOptions(c_atts)

    # Streamline --------------------------------
    sl_atts = StreamlineAttributes()
    sl_atts.sourceType = sl_atts.SpecifiedBox # center the values at cell corners
    sl_atts.stepLength = sline_StepLength
    sl_atts.maxTime = sline_MaxSteps
    sl_atts.boxExtents = (sline_BoxX[0], sline_BoxX[1],
                          kf_SlineBoxYLo[0], kf_SlineBoxYHi[0],
                          kf_SlineBoxZLo[0], kf_SlineBoxZHi[0])
    sl_atts.pointDensity = sline_PointDensity
    sl_atts.displayMethod = sl_atts.Tubes
    sl_atts.showStart = sline_ShowStart
    sl_atts.radius = kf_SlineRadius[0]
    sl_atts.coloringMethod = sl_atts.ColorBySpeed
    SetPlotOptions(sl_atts)

    # Set up operator(s)
    if num_dims == "3d":

        SetActivePlots((0, 1)) # Make both plots active for Transform op (seems needed here to make both active)
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
        SetOperatorOptions(xform_atts)

        SetActivePlots(0) # Make only Contour plot active for Slice op
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)

        SetActivePlots(1) # Make only Streamline plot active for Resample op
        AddOperator("Resample")
        res_atts = ResampleAttributes()
        res_atts.defaultVal = 0.0
        res_atts.width = resample_SamplesXYZ[0]
        res_atts.height = resample_SamplesXYZ[1]
        res_atts.depth = resample_SamplesXYZ[2]
        res_atts.minX = resample_StartEndX[0]
        res_atts.maxX = resample_StartEndX[1]
        res_atts.minY = kf_ResampleStartY[0]
        res_atts.maxY = kf_ResampleEndY[0]
        res_atts.minZ = kf_ResampleStartZ[0]
        res_atts.maxZ = kf_ResampleEndZ[0]
        SetOperatorOptions(res_atts)
        AddOperator("Clip")
        clip_atts = ClipAttributes()
        clip_atts.funcType = clip_atts.Plane
        clip_atts.plane1Normal = clip_Plane1Normal
        clip_atts.plane1Origin = clip_Plane1Origin
        clip_atts.plane1Status = 1
        clip_atts.plane2Status = 0
        clip_atts.plane3Status = 0
        clip_atts.planeInverse = 0
        SetOperatorOptions(clip_atts)

        SetActivePlots((0, 1)) # Make both plots active

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window

    return sl_atts


#================================================================================
# Set up plots and operators for contour & streamlines
def set_plots_ops_7a():
    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Contour", visit_vname)
    AddPlot("Streamline", "velocity_3d")
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
    l_atts.direction = (-1, 0, -1)
    SetLight(0, l_atts)

    # Contour --------------------------------
    c_atts = ContourAttributes()
    c_atts.contourMethod = c_atts.Value # center the values at cell corners
    c_atts.contourValue = 0.5
    c_atts.lineWidth = 3
#    c_atts.limitsMode = 1 # limits a function of 0) original data or 1) current plot
#    if commonrange == "True":
#        c_atts.min = highmin # lower bound of range that's colormapped
#        c_atts.max = lowmax # upper bound of range that's colormapped
#        c_atts.minFlag = 1 # turn on
#        c_atts.maxFlag = 1 # turn on
    SetPlotOptions(c_atts)

    # Streamline --------------------------------
    sl_atts = StreamlineAttributes()
    sl_atts.sourceType = sl_atts.SpecifiedBox # center the values at cell corners
    sl_atts.stepLength = sline_StepLength
    sl_atts.maxTime = sline_MaxSteps
    sl_atts.boxExtents = (sline_BoxX[0], sline_BoxX[1],
                          kf_SlineBoxYLo[0], kf_SlineBoxYHi[0],
                          kf_SlineBoxZLo[0], kf_SlineBoxZHi[0])
    sl_atts.pointDensity = sline_PointDensity
    sl_atts.displayMethod = sl_atts.Tubes
    sl_atts.showStart = sline_ShowStart
    sl_atts.radius = kf_SlineRadius[0] * radius_factor
    sl_atts.radius = kf_SlineRadius[0]
    sl_atts.coloringMethod = sl_atts.ColorBySpeed
    SetPlotOptions(sl_atts)

    # Set up operator(s)
    if num_dims == "3d":

        SetActivePlots((0, 1)) # Make both plots active for Transform op (seems needed here to make both active)
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
        SetOperatorOptions(xform_atts)

        SetActivePlots(0) # Make only Contour plot active for Slice op
        AddOperator("Slice")
        s_atts = SliceAttributes()
        s_atts.project2d = 0
        s_atts.normal = (1,0,0)
        s_atts.upAxis = (0,1,0)
        SetOperatorOptions(s_atts)

        SetActivePlots(1) # Make only Streamline plot active for Resample op
        AddOperator("Clip")
        clip_atts = ClipAttributes()
        clip_atts.funcType = clip_atts.Plane
        clip_atts.plane1Normal = clip_Plane1Normal
        clip_atts.plane1Origin = clip_Plane1Origin
        clip_atts.plane1Status = 1
        clip_atts.plane2Status = 0
        clip_atts.plane3Status = 0
        clip_atts.planeInverse = 0
        SetOperatorOptions(clip_atts)

        SetActivePlots((0, 1)) # Make both plots active

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window

    return sl_atts


#================================================================================
# Set up plots and operators for contour & streamlines
def set_plots_ops_7b():
    DisableRedraw() # disable redrawing for each of many operations

    # Set up Resample operator
    if num_dims == "3d":

        SetActivePlots(1) # Make only Streamline plot active for Resample op
        AddOperator("Resample")
        res_atts = ResamplePluginAttributes()
        res_atts.defaultValue = 0.0
        res_atts.tieResolver = res_atts.random
        res_atts.is3D = True
        res_atts.samplesX = resample_SamplesXYZ[0]
        res_atts.samplesY = resample_SamplesXYZ[1]
        res_atts.samplesZ = resample_SamplesXYZ[2]
        res_atts.startX = resample_StartEndX[0]
        res_atts.endX = resample_StartEndX[1]
        res_atts.startY = kf_ResampleStartY[0]
        res_atts.endY = kf_ResampleEndY[0]
        res_atts.startZ = kf_ResampleStartZ[0]
        res_atts.endZ = kf_ResampleEndZ[0]
        SetOperatorOptions(res_atts)

        SetActivePlots((0, 1)) # Make both plots active

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window

    return res_atts


#================================================================================
# Set up plots and operators for contour & streamlines
def set_plots_8():
    # ---------------------------------------------------------------------------
    # Set up plot(s)
    visit_vname = "mesh_blockandlevel/" + vname
    AddPlot("Contour", visit_vname) # For isosurface
    AddPlot("Contour", visit_vname) # For edge of clipped isosurface
#    AddPlot("Pseudocolor", visit_vname) # For pseudocolor of interior on faces exposed by clipping
    AddPlot("Vector", "velocity_3d") # Hedgehogs
    DisableRedraw() # disable redrawing for each of many operations

    l_atts = LightAttributes()
    l_atts.direction = (-1, -1, -1)
    SetLight(0, l_atts)

    # Contour surface --------------------------------
    SetActivePlots(0) # Isosurface plot
    c0_atts = ContourAttributes()
    c0_atts.contourMethod = c0_atts.Value # center the values at cell corners
    c0_atts.contourValue = 0.5
    c0_atts.lineWidth = 3
    SetPlotOptions(c0_atts)
    SetActivePlots((0, 1, 2)) # All plots
#    HideActivePlots()

    # Contour surface edge --------------------------------
    SetActivePlots(1) # Edge plot
    c1_atts = ContourAttributes()
    c1_atts.contourMethod = c1_atts.Value # center the values at cell corners
    c1_atts.contourValue = 0.5
#    c1_atts.contourValue = 0.25
    c1_atts.lineWidth = 3
    c1_atts.colorType = c1_atts.ColorBySingleColor
    c1_atts.singleColor = (255, 255, 0, 255)
#    c1_atts.wireframe = 1
    SetPlotOptions(c1_atts)
    SetActivePlots((0, 1, 2)) # All plots

#     # Pseudocolor ---------------------------
#     SetActivePlots(2) # Pseudocolor plot
#     pc_atts = PseudocolorAttributes()
#     pc_atts.centering = pc_atts.Nodal # center the values at cell corners
#     pc_atts.opacity = 1.0
# #    pc_atts.limitsMode = 1 # limits a function of 0) original data or 1) current plot
#     if commonrange == "True":
#         pc_atts.min = highmin # lower bound of range that's colormapped
#         pc_atts.max = lowmax # upper bound of range that's colormapped
#         pc_atts.minFlag = 1 # turn on
#         pc_atts.maxFlag = 1 # turn on
#     SetPlotOptions(pc_atts)
#     SetActivePlots((0, 1, 2)) # All plots

    # Vector ---------------------------------
    SetActivePlots(2) # Vector plot
    vec_atts = VectorAttributes()
#    vec_atts.useStride = 0
#    numvecs = numvecs_start
#    vec_atts.nVectors = 70
    vec_atts.useStride = 1
    vec_atts.stride = int(kf_vecstride[0])
    vec_atts.scale = kf_hhogscale[0]
    vec_atts.scaleByMagnitude = 0
    vec_atts.headSize = 0.22
    vec_atts.lineStem = 0
    vec_atts.highQuality = 1
    vec_atts.stemWidth = 0.19
    SetPlotOptions(vec_atts)
    SetActivePlots((0, 1, 2)) # All plots

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window

    return vec_atts

#================================================================================
# Set up plots and operators for contour & streamlines
def set_ops_8():
    DisableRedraw() # disable redrawing for each of many operations

    # Set up Resample operator
    if num_dims == "3d":

        SetActivePlots((0, 1)) # All but Vector plot
#        SetActivePlots(1) # Make only Streamline plot active for Resample op
        AddOperator("Clip")
        clip_atts = ClipAttributes()
        clip_atts.funcType = clip_atts.Plane
        clip_atts.plane1Normal = clip1_Plane1Normal
        clip_atts.plane1Origin = clip1_Plane1Origin
        clip_atts.plane2Normal = clip1_Plane2Normal
        clip_atts.plane2Origin = clip1_Plane2Origin
        clip_atts.plane1Status = 1
        clip_atts.plane2Status = 1
        clip_atts.plane3Status = 0
        clip_atts.planeInverse = 0
        SetOperatorOptions(clip_atts)
        SetActivePlots((0, 1, 2)) # All plots

        SetActivePlots((0, 1, 2)) # All plots
        AddOperator("Transform")
        xform_atts = TransformAttributes()
        xform_atts.rotateAxis = (1, 0, 0)
        xform_atts.rotateAmount = -90
        xform_atts.doRotate = 1
        SetOperatorOptions(xform_atts)
        SetActivePlots((0, 1, 2)) # All plots

        SetActivePlots(1) # Make 2nd Contour plot active
        AddOperator("ExternalSurface")
        es_atts = ExternalSurfaceAttributes()
        SetOperatorOptions(es_atts)
        SetActivePlots((0, 1, 2)) # All plots

        SetActivePlots(2) # Make Vector plot active

        AddOperator("ThreeSlice")
        ts_atts = ThreeSliceAttributes()
        SetOperatorOptions(ts_atts)
        AddOperator("Clip")
        clip2_atts = ClipAttributes()
        clip2_atts.funcType = clip2_atts.Plane
        clip2_atts.plane1Normal = clip2_Plane1Normal
        clip2_atts.plane1Origin = clip2_Plane1Origin
        clip2_atts.plane2Normal = clip2_Plane2Normal
        clip2_atts.plane2Origin = clip2_Plane2Origin
        clip2_atts.plane1Status = 1
        clip2_atts.plane2Status = 0 # 1
        clip2_atts.plane3Status = 0
        clip2_atts.planeInverse = 1 # 1
        SetOperatorOptions(clip2_atts)
#         AddOperator("Slice")
#         s_atts = SliceAttributes()
#         SetOperatorOptions(s_atts)

        AddOperator("Threshold")
        t_atts = ThresholdAttributes()
#        t_atts.defaultVarName = "mesh_blockandlevel/rpv1"
        t_atts.listedVarNames = "mesh_blockandlevel/rpv1"
        t_atts.lowerBounds = 0.45
        print "[set_ops_8()] - t_atts.listedVarNames: ", t_atts.listedVarNames
#        t_atts.upperBounds = 0.9999
        SetOperatorOptions(t_atts)
        SetActivePlots((0, 1, 2)) # All plots

    DrawPlots()
    RedrawWindow() # (re)enable redrawing of vis window


#================================================================================
# Select one of the refinement levels for display
def set_ref_level():
    # ---------------------------------------------------------------------------
    # Select refinement level
    silr = SILRestriction()
    cats = silr.Categories()
    l_cats = len(cats)
#    print "Length of cats: ", l_cats
#    print "silr.Categories: ", cats
#    print "silr.NumCategories: ", silr.NumCategories()
#    print "silr.NumSets: ", silr.NumSets()

    sic = silr.SetsInCategory(cats[1])
#    print "sic: ", sic
    l_sic = len(sic)
#    print "Length of sic: ", l_sic
#    print "Category: ", cats[1]
#    print "Sets in Category: ", sic

    silr.TurnOffAll()
    SetPlotSILRestriction(silr)

    silr.TurnOnSet(sic[l_sic-1]) # silr.TurnOnSet(sic[l_sic-2]), silr.TurnOnSet(sic[l_sic-3])
#    silr.TurnOnSet(sic[(l_sic-1)/2])
#    silr.TurnOnSet( sic[ ((l_sic-1)/2) + 1 ] )
#    silr.TurnOnSet( sic[ 7 ] )
    print "----------------------------------------------------------------------"
    print "[set_ref_level()] - Isolating data to 1 refinement level: ", sic[ l_sic-1 ] # sic[ l_sic-3 ]
#    print "[set_ref_level()] - Isolating data to 1 refinement level: ", sic[ ((l_sic-1)/2) + 1 ]
    print "----------------------------------------------------------------------"
#    silr.TurnOnSet(sic[2])
#    print "(l_sic-1)/2: ", (l_sic-1)/2
#    print "l_sic-2: ", l_sic-2

    SetPlotSILRestriction(silr)
#    print "Turned off all and turned on only this subset: ", sic[l_sic-1]
#    print "Turned off all and turned on only this subset: ", sic[l_sic-2]
#    print "Turned off all and turned on only this subset: ", sic[(l_sic-1)/2]




#================================================================================
# One of several rendering loops that can be called by snap()
def snaploop_1():

    slider = CreateAnnotationObject("TimeSlider")

    for state in range(1, TimeSliderGetNStates(), state_incr):  #  Skip first file
        SetTimeSliderState(state)
        Query("Time")
        newtime = GetQueryOutputValue()*1.0
        slider.text = "Time = %g seconds" % newtime
        SaveWindow()

#        SetActivePlots(0) # Make only Contour plot active for Query
#        Query("Centroid") # get centroid of Contour plot
#        centroid = GetQueryOutputValue()
#        SetActivePlots((0, 1)) # Make both plots active

        #------------------------------------------------------------------------
#        zoom += zoom_incr
#        view3d.imageZoom = zoom
#        focus = (
#            focus[0] + focus_incr[0],
#            focus[1] + focus_incr[1],
#            focus[2] + focus_incr[2]
#            )
#        view3d.focus = focus

        #------------------------------------------------------------------------
#        pan = (
#            pan[0] + pan_incr[0],
#            pan[1] + pan_incr[1]
#            )
#        view3d.imagePan = pan

#        print "-----------------------------------------------------------------------"
#        print "[snap()] - zoom: ", zoom
#        print "[snap()] - pan: ", pan
#        print "-----------------------------------------------------------------------"

        #------------------------------------------------------------------------
#        SetView3D(view3d)
#        set_viewscope()

        #------------------------------------------------------------------------
#        numvecs += numvecs_incr
#        vecscale += vecscale_incr
#        vec_atts.nVectors = numvecs
#        vec_atts.scale = vecscale
#        SetPlotOptions(vec_atts)

        #------------------------------------------------------------------------
#        streamplaneradius += streamplaneradius_incr
#        sl_atts.planeRadius = streamplaneradius
##        streamplaneorigin = (
##            streamplaneorigin[0] + streamplaneorigin_incr[0],
##            streamplaneorigin[1] + streamplaneorigin_incr[1],
##            streamplaneorigin[2] + streamplaneorigin_incr[2]
##            )
##        sl_atts.planeOrigin = (streamplaneorigin[0], streamplaneorigin[1], streamplaneorigin[2])
#        sl_atts.planeOrigin = (centroid[0], centroid[1], centroid[2])
##        sl_atts.planeOrigin = streamplaneorigin
#        SetPlotOptions(sl_atts)


#================================================================================
# One of several rendering loops that can be called by snap()
def snaploop_2(vec_atts):

#    print "-----------------------------------------------------------------------"
#    print "[snaploop_2()] - vec_atts: ", vec_atts
#    print "-----------------------------------------------------------------------"

## Create a tuple of camera values and x values. The x values determine where in [0,1] the control points occur.
#cpts = (c0, c1, c2, c3, c4, c5, c6)
#x=[]
#for i in range(7):
#    x = x + [float(i) / float(6.)]
## Animate the view using EvalCubicSpline.
#nsteps = 100
#for i in range(nsteps):
#    t = float(i) / float(nsteps - 1)
#    c = EvalCubicSpline(t, x, cpts)
#    c.nearPlane = -34.461
#    c.farPlane = 34.461
#    SetView3D(c)

    weights = []
#     len_ts = len(kf_timesteps)
#     range_kfts = kf_timesteps[-1] - kf_timesteps[0]
#     for i in range(len_ts):
#         dist = kf_timesteps[i] - kf_timesteps[0] # distance into kf_timesteps
#         weights = weights + [dist / range_kfts]
    range_kfts = kf_timesteps[-1] - kf_timesteps[0]
    for kf_step in kf_timesteps:
        dist = kf_step - kf_timesteps[0] # distance into kf_timesteps
        weights = weights + [dist / range_kfts]
    print "-----------------------------------------------------------------------"
    print "[snaploop_2()] - weights: ", weights
    print "-----------------------------------------------------------------------"

    slider = CreateAnnotationObject("TimeSlider")

    if len(timesteps) != TimeSliderGetNStates():
        print "-------------------------------------------------------------------"
        print "[snaploop_2()] - Number of timesteps NE number of slider states"
        print "[snaploop_2()] -    Exiting..."
        print "-------------------------------------------------------------------"
        sys.exit()

    range_ts = timesteps[-1] - timesteps[0]

#    for state in range(1, TimeSliderGetNStates(), state_incr):  #  Skip first file
    for state in range(0, TimeSliderGetNStates(), state_incr):  #  DO NOT skip first file

        print "-------------------------------------------------------------------"
        print "[snaploop_2()] - Current file: ", fnames[state]
        print "-------------------------------------------------------------------"

        t = (timesteps[state] - timesteps[0]) / range_ts
        this_vecstride = EvalCubicSpline(t, weights, kf_vecstride)
        print "-------------------------------------------------------------------"
        print "[snaploop_2()] - This vector stride: ", this_vecstride
        print "-------------------------------------------------------------------"
        this_hhogscale = EvalCubicSpline(t, weights, kf_hhogscale)
        print "-------------------------------------------------------------------"
        print "[snaploop_2()] - This hedgehog scale: ", this_hhogscale
        print "-------------------------------------------------------------------"
        vec_atts.stride = int(this_vecstride)
        vec_atts.scale = this_hhogscale
#        SetPlotOptions(vec_atts)

        SetTimeSliderState(state)
        Query("Time")
        newtime = GetQueryOutputValue()*1.0
        slider.text = "Time = %g seconds" % newtime
        SaveWindow()









#        SetActivePlots(0) # Make only Contour plot active for Query
#        Query("Centroid") # get centroid of Contour plot
#        centroid = GetQueryOutputValue()
#        SetActivePlots((0, 1)) # Make both plots active

        #------------------------------------------------------------------------
#        zoom += zoom_incr
#        view3d.imageZoom = zoom
#        focus = (
#            focus[0] + focus_incr[0],
#            focus[1] + focus_incr[1],
#            focus[2] + focus_incr[2]
#            )
#        view3d.focus = focus

        #------------------------------------------------------------------------
#        pan = (
#            pan[0] + pan_incr[0],
#            pan[1] + pan_incr[1]
#            )
#        view3d.imagePan = pan

#        print "-----------------------------------------------------------------------"
#        print "[snap()] - zoom: ", zoom
#        print "[snap()] - pan: ", pan
#        print "-----------------------------------------------------------------------"

        #------------------------------------------------------------------------
#        SetView3D(view3d)
#        set_viewscope()

        #------------------------------------------------------------------------
#        numvecs += numvecs_incr
#        vecscale += vecscale_incr
#        vec_atts.nVectors = numvecs
#        vec_atts.scale = vecscale
#        SetPlotOptions(vec_atts)

        #------------------------------------------------------------------------
#        streamplaneradius += streamplaneradius_incr
#        sl_atts.planeRadius = streamplaneradius
##        streamplaneorigin = (
##            streamplaneorigin[0] + streamplaneorigin_incr[0],
##            streamplaneorigin[1] + streamplaneorigin_incr[1],
##            streamplaneorigin[2] + streamplaneorigin_incr[2]
##            )
##        sl_atts.planeOrigin = (streamplaneorigin[0], streamplaneorigin[1], streamplaneorigin[2])
#        sl_atts.planeOrigin = (centroid[0], centroid[1], centroid[2])
##        sl_atts.planeOrigin = streamplaneorigin
#        SetPlotOptions(sl_atts)


#================================================================================
# One of several rendering loops that can be called by snap()
def snaploop_3(sl_atts, res_atts):
    weights = []
    range_kfts = kf_timesteps[-1] - kf_timesteps[0]
    for kf_step in kf_timesteps:
        dist = kf_step - kf_timesteps[0] # distance into kf_timesteps
        weights = weights + [dist / range_kfts]
    print "-----------------------------------------------------------------------"
    print "[snaploop_3()] - weights: ", weights
    print "-----------------------------------------------------------------------"

    slider = CreateAnnotationObject("TimeSlider")

    if len(timesteps) != TimeSliderGetNStates():
        print "-------------------------------------------------------------------"
        print "[snaploop_3()] - Number of timesteps NE number of slider states"
        print "[snaploop_3()] -    Exiting..."
        print "-------------------------------------------------------------------"
        sys.exit()

    range_ts = timesteps[-1] - timesteps[0]

#    for state in range(1, TimeSliderGetNStates(), state_incr):  #  Skip first file
    for state in range(0, TimeSliderGetNStates(), state_incr):  #  DO NOT skip first file
        print "-------------------------------------------------------------------"
        print "[snaploop_3()] - Current file: ", fnames[state]
        print "-------------------------------------------------------------------"

        t = (timesteps[state] - timesteps[0]) / range_ts
        this_SlineRadius = EvalCubicSpline(t, weights, kf_SlineRadius)
        this_SlineBoxYLo = EvalCubicSpline(t, weights, kf_SlineBoxYLo)
        this_SlineBoxYHi = EvalCubicSpline(t, weights, kf_SlineBoxYHi)
        this_SlineBoxZLo = EvalCubicSpline(t, weights, kf_SlineBoxZLo)
        this_SlineBoxZHi = EvalCubicSpline(t, weights, kf_SlineBoxZHi)
        this_ResampleStartY = EvalCubicSpline(t, weights, kf_ResampleStartY)
        this_ResampleEndY = EvalCubicSpline(t, weights, kf_ResampleEndY)
        this_ResampleStartZ = EvalCubicSpline(t, weights, kf_ResampleStartZ)
        this_ResampleEndZ = EvalCubicSpline(t, weights, kf_ResampleEndZ)

        sl_atts.boxExtents = (sline_BoxX[0], sline_BoxX[1],
                              this_SlineBoxYLo, this_SlineBoxYHi,
                              this_SlineBoxZLo, this_SlineBoxZHi)
        sl_atts.radius = this_SlineRadius * radius_factor
#        sl_atts.radius = this_SlineRadius
        SetPlotOptions(sl_atts)

        res_atts.startY = this_ResampleStartY
        res_atts.endY = this_ResampleEndY
        res_atts.startZ = this_ResampleStartZ
        res_atts.endZ = this_ResampleEndZ
        SetOperatorOptions(res_atts)

        SetTimeSliderState(state)
        Query("Time")
        newtime = GetQueryOutputValue()*1.0
        slider.text = "Time = %g seconds" % newtime
        SaveWindow()


#================================================================================
# One of several rendering loops that can be called by snap()
#def snaploop_8():
def snaploop_8(vec_atts):
    weights = []
    range_kfts = kf_timesteps[-1] - kf_timesteps[0]
    for kf_step in kf_timesteps:
        dist = kf_step - kf_timesteps[0] # distance into kf_timesteps
        weights = weights + [dist / range_kfts]

    slider = CreateAnnotationObject("TimeSlider")

    if len(timesteps) != TimeSliderGetNStates():
        print "-------------------------------------------------------------------"
        print "[snaploop_8()] - Number of timesteps NE number of slider states"
        print "[snaploop_8()] -    Exiting..."
        print "-------------------------------------------------------------------"
        sys.exit()

    range_ts = timesteps[-1] - timesteps[0]

    for state in range(0, TimeSliderGetNStates(), state_incr):  #  DO NOT skip first file
        if timesteps[state] > 165.0:
            print "-------------------------------------------------------------------"
            print "[snaploop_8()] - Current file: ", fnames[state]
            print "-------------------------------------------------------------------"

            t = (timesteps[state] - timesteps[0]) / range_ts
            this_vecstride = EvalCubicSpline(t, weights, kf_vecstride)
            print "-------------------------------------------------------------------"
            print "[snaploop_2()] - This vector stride: ", this_vecstride
            print "-------------------------------------------------------------------"
            this_hhogscale = EvalCubicSpline(t, weights, kf_hhogscale)
            print "-------------------------------------------------------------------"
            print "[snaploop_2()] - This hedgehog scale: ", this_hhogscale
            print "-------------------------------------------------------------------"
            vec_atts.stride = int(this_vecstride)
            vec_atts.scale = this_hhogscale
            SetPlotOptions(vec_atts)

            SetTimeSliderState(state)
            Query("Time")
            newtime = GetQueryOutputValue()*1.0
            slider.text = "Time = %g seconds" % newtime
            SaveWindow()


#================================================================================
# Snap images
def snap():

    global zoom, zoom_start, zoom_stop, zoom_incr
    global pan, pan_start, pan_stop, pan_incr
    global state_incr
    global numvecs, numvecs_start, vecscale, vecscale_start
    global streamplaneorigin, streamplaneorigin_start, streamplaneorigin_stop, streamplaneorigin_incr
    global streamplaneradius, streamplaneradius_start, streamplaneradius_stop, streamplaneradius_incr

    # ---------------------------------------------------------------------------
    # Find intersection of ranges of variable values over time series
    if commonrange == "True":
        print "[snap()] - Calling range_intersection()"
        range_intersection()
        print "[snap()] - All minima: ", minima
        print "[snap()] - All maxima: ", maxima
        print "[snap()] - Lowest maximum: ", lowmax
        print "[snap()] - Highest minimum: ", highmin

    # ---------------------------------------------------------------------------
#    if engine_loc == "davinci":
    if dbase_loc == "davinci":
        # Launch parallel compute engine remotely on davinci
        print "------------------------------------------------------------------"
        print "[snap()] - Launching parallel compute engine on davinci"
        print "------------------------------------------------------------------"
        args = ("-np", "12") # ("-np", "8", "-nn", "8")
        OpenComputeEngine("davinci.nersc.gov", args)
#    elif engine_loc == "local":
    elif dbase_loc == "local":
        print "------------------------------------------------------------------"
        print "[snap()] - Compute engine will run locally"
        print "------------------------------------------------------------------"

    # ---------------------------------------------------------------------------
    # Open time series of files
##    if engine_loc == "davinci":
#    if dbase_loc == "davinci" and client_loc == "local":
#        dbasename = rmt_machine_pfx + pf_path + pfn_prefix + "*" + " database"
##    elif engine_loc == "local":
#    elif dbase_loc == "local":
#        dbasename = pf_path + pfn_prefix + "*" + " database"

    dbasename = pf_path + pfn_prefix + "*" + " database"
    if dbase_loc == "davinci" and client_loc == "local":
        dbasename = rmt_machine_pfx + dbasename

    OpenDatabase(dbasename)
    print "Number of time states: %d" % GetDatabaseNStates()

#    set_viewscope() # Looks like this ought to be called AFTER SetView3D()
  
    # ---------------------------------------------------------------------------
    # Set up plots and operators
#    set_plots_ops_1() # pseudocolor & vectors
#    set_plots_ops_2() # contour
#    set_plots_ops_3() # vectors
#    set_plots_ops_4() # streamlines
#    sl_atts = set_plots_ops_5() # streamlines & contour slice
#    vec_atts = set_plots_ops_6() # hedgehogs & contour slice

    if movie_id == "sline_slice":
        sl_atts = set_plots_ops_7a() # streamlines & contour slice
        res_atts = set_plots_ops_7b() # streamlines & contour slice
    elif movie_id =="hhogs_isosurf_clipped":
#        set_plots_8() # Two Contours & Vector
        vec_atts = set_plots_8() # Two Contours & Vector
        set_ops_8() # One Clip of both Contours; Thresholded ThreeSlice of Vector
    print "----------------------------------------------------------------------"
    print "[snap()] - Just set up plots (and, maybe, operators)"
    print "----------------------------------------------------------------------"

#    ToggleMaintainViewMode()
#    ToggleMaintainDataMode()

    # ---------------------------------------------------------------------------
    # Select refinement level to display
    set_ref_level()
    print "----------------------------------------------------------------------"
    print "[snap()] - Just set ref level"
    print "----------------------------------------------------------------------"

#    an_atts = GetAnnotationAttributes()
#    an_atts.xGridLines = 1
#    an_atts.yGridLines = 1
#    an_atts.zGridLines = 1
#    an_atts.backgroundColor = (70, 83, 235) # (65, 83, 230)
#    an_atts.foregroundColor = (255, 255, 255)
#    SetAnnotationAttributes(an_atts)

    # ---------------------------------------------------------------------------
    # Set view parameters
    if num_dims == "3d":
        view3d = GetView3D()
        view3d.viewNormal = camvec
        view3d.viewUp = (tuple(upvec))

#        zoom = zoom_start
#        view3d.imageZoom = zoom
#        focus = focus_start
#        view3d.focus = focus
#        pan = pan_start
#        view3d.imagePan = pan

        SetView3D(view3d)
    elif num_dims == "2d":
        view2d = View2DAttributes()
        print "[snap()] - view port coords: ", view2d.viewportCoords
        print "[snap()] - window coords: ", view2d.windowCoords
        view2d.viewportCoords = (0.4, 0.95, 0.15, 0.95)
        view2d.windowCoords = (0.0, 2.46e9, -2.46e9, 2.46e9)
        SetView2D(view2d)

    set_viewscope()

    # ---------------------------------------------------------------------------
    # Set saved-image parameters
    swa = SaveWindowAttributes()
    swa.format = swa.PNG
#    imagefilename = if_path + ifn_prefix + "%04d"%(i,)
#    swa.fileName = imagefilename
    swa.family = 1 # saving series ("family") of images
    SetSaveWindowAttributes(swa)

#    zoom_incr = (zoom_stop - zoom_start)/(TimeSliderGetNStates() - 1)
#    focus_incr = (
#        state_incr * ((zoom_stop - zoom_start)/((TimeSliderGetNStates() - 1) - 1))
#        (focus_stop[0] - focus_start[0])/(TimeSliderGetNStates() - 1),
#        (focus_stop[1] - focus_start[1])/(TimeSliderGetNStates() - 1),
#        (focus_stop[2] - focus_start[2])/(TimeSliderGetNStates() - 1)
#        )
#    print "[snap()] - TimeSliderGetNStates(): ", TimeSliderGetNStates()
    zoom_incr = state_incr * ((zoom_stop - zoom_start)/((TimeSliderGetNStates() - 1) - 1))
#    focus_incr = (
#        state_incr * ((focus_stop[0] - focus_start[0])/((TimeSliderGetNStates() - 1) - 1)),
#        state_incr * ((focus_stop[1] - focus_start[1])/((TimeSliderGetNStates() - 1) - 1)),
#        state_incr * ((focus_stop[2] - focus_start[2])/((TimeSliderGetNStates() - 1) - 1))
#        )
    pan_incr = (
        state_incr * ((pan_stop[0] - pan_start[0])/((TimeSliderGetNStates() - 1) - 1)),
        state_incr * ((pan_stop[1] - pan_start[1])/((TimeSliderGetNStates() - 1) - 1))
        )
#    print "[snap()] - zoom_incr: ", zoom_incr
##    print "[snap()] - focus_incr[0] = %d * ((%f - %f)/((%d - 1) - 1))" % (state_incr, focus_stop[0], focus_start[0], TimeSliderGetNStates())
##    print "[snap()] - focus_incr: ", focus_incr
#    print "[snap()] - pan_incr: ", pan_incr

#    numvecs_incr = state_incr * ((numvecs_stop - numvecs_start)/((TimeSliderGetNStates() - 1) - 1))
    numvecs_incr = ((state_incr * (numvecs_stop - numvecs_start))/((TimeSliderGetNStates() - 1) - 1))
    vecscale_incr = state_incr * ((vecscale_stop - vecscale_start)/((TimeSliderGetNStates() - 1) - 1))
#    print "-----------------------------------------------------------------------"
#    print "[snap()] - numvecs_incr: ", numvecs_incr
#    print "[snap()] - vecscale_incr: ", vecscale_incr
#    print "-----------------------------------------------------------------------"

    streamplaneorigin_incr = (
        state_incr * ((streamplaneorigin_stop[0] - streamplaneorigin_start[0])/((TimeSliderGetNStates() - 1) - 1)),
        state_incr * ((streamplaneorigin_stop[1] - streamplaneorigin_start[1])/((TimeSliderGetNStates() - 1) - 1)),
        state_incr * ((streamplaneorigin_stop[2] - streamplaneorigin_start[2])/((TimeSliderGetNStates() - 1) - 1))
        )
    streamplaneradius_incr = state_incr * ((streamplaneradius_stop - streamplaneradius_start)/((TimeSliderGetNStates() - 1) - 1))
    print "-----------------------------------------------------------------------"
    print "[snap()] - streamplaneorigin_incr: ", streamplaneorigin_incr
    print "[snap()] - streamplaneradius_incr: ", streamplaneradius_incr
    print "-----------------------------------------------------------------------"

#    snaploop_1()
#    snaploop_2(vec_atts)
    if movie_id == "sline_slice":
        snaploop_3(sl_atts, res_atts)
    elif movie_id =="hhogs_isosurf_clipped":
#        snaploop_8()
        snaploop_8(vec_atts)

    sys.exit()

#---------\---------\---------\---------\---------\---------\---------\---------\
#---------/---------/---------/---------/---------/---------/---------/---------/
#    # ---------------------------------------------------------------------------
#    # 1) Define time slider
#    # 2) Loop over time steps and
#    #    a) change time
#    #    b) stuff new time into slider caption
#    #    c) snap image
##    for state in range(TimeSliderGetNStates()):
#    slider = CreateAnnotationObject("TimeSlider")
#    for state in range(1, TimeSliderGetNStates(), state_incr):  #  Skip first file
#        print "------------------------------------------------------------------"
#        print "[snap()]- loop - vec_atts: ", vec_atts
#        print "------------------------------------------------------------------"
#        SetTimeSliderState(state)
#        Query("Time")
##        slider = CreateAnnotationObject("TimeSlider")
##        print "Type of slider: ", type(slider)
#        newtime = GetQueryOutputValue()*1.0
##        print "Type of str(newtime): ", type(str(newtime))
#        slider.text = "Time = %g seconds" % newtime
#        SaveWindow()
#
#        SetActivePlots(0) # Make only Contour plot active for Query
#        Query("Centroid") # get centroid of Contour plot
#        print "=================================================================="
#        print "[snap()]- Just queried Centroid of Contour"
#        print GetQueryOutputString()
#        print "[snap()]- Just printed query string"
#        print GetQueryOutputValue()
#        centroid = GetQueryOutputValue()
#        print "[snap()]- Just printed query value"
#        print "=================================================================="
#        SetActivePlots((0, 1)) # Make both plots active
#
#        #------------------------------------------------------------------------
##        zoom += zoom_incr
##        view3d.imageZoom = zoom
##        focus = (
##            focus[0] + focus_incr[0],
##            focus[1] + focus_incr[1],
##            focus[2] + focus_incr[2]
##            )
##        view3d.focus = focus
#
#        #------------------------------------------------------------------------
##        pan = (
##            pan[0] + pan_incr[0],
##            pan[1] + pan_incr[1]
##            )
##        view3d.imagePan = pan
#
##        print "-----------------------------------------------------------------------"
##        print "[snap()] - zoom: ", zoom
##        print "[snap()] - pan: ", pan
##        print "-----------------------------------------------------------------------"
#
#        #------------------------------------------------------------------------
##        SetView3D(view3d)
##        set_viewscope()
#
#        #------------------------------------------------------------------------
##        numvecs += numvecs_incr
##        vecscale += vecscale_incr
##        vec_atts.nVectors = numvecs
##        vec_atts.scale = vecscale
##        SetPlotOptions(vec_atts)
#
#        #------------------------------------------------------------------------
##        streamplaneradius += streamplaneradius_incr
##        sl_atts.planeRadius = streamplaneradius
###        streamplaneorigin = (
###            streamplaneorigin[0] + streamplaneorigin_incr[0],
###            streamplaneorigin[1] + streamplaneorigin_incr[1],
###            streamplaneorigin[2] + streamplaneorigin_incr[2]
###            )
###        sl_atts.planeOrigin = (streamplaneorigin[0], streamplaneorigin[1], streamplaneorigin[2])
##        sl_atts.planeOrigin = (centroid[0], centroid[1], centroid[2])
###        sl_atts.planeOrigin = streamplaneorigin
##        SetPlotOptions(sl_atts)
#---------\---------\---------\---------\---------\---------\---------\---------\
#---------/---------/---------/---------/---------/---------/---------/---------/


#--------------------------------------------------------------------------------
# This is the script's "main()"
# Make sure the "up" and "normal" vectors of view plane don not coincide
if camvec == tuple(upvec_plus_y):
    upvec = tuple(upvec_minus_z)
    print "[main] - camvec == upvec_plus_y, so upvec will NOT be default: ", upvec
else:
    upvec = tuple(upvec_plus_y)
    print "[main] - camvec != upvec_plus_y, so upvec WILL be default: ", upvec
snap()
sys.exit()



#    # ---------------------------------------------------------------------------
#    # TEST STRING ARITHMETIC
#    tstr = "middle"
#    print "string start", tstr
#    tstr = "prefix" + tstr + "suffix"
#    print "string after arithmetic", tstr
#    sys.exit()
#    # ---------------------------------------------------------------------------
