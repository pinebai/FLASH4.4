#!/usr/bin/env python

# CHANGELOG:
#  2016-10-18 Updated "print" calls for python 3.x and 2.x cross-compatibility. Replaced "xrange" with "range" -Scott Feister
#             Updated getWhereList for cross-compatibility

# import visit_writer

import tables
import numpy as np
from sys import argv

import warnings
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

# Load the data:
def write_vtk(filename):
    print("Processing file:" + filename)

    h5file = tables.openFile(filename, "r")

    # Read nstep and time from plot file:
    n = "%-80s" % "nstep"
    t = h5file.getNode("/integer scalars")
    if int(tables.__version__.split(".")[0]) >= 3: # Pytables 3.x syntax
        wl = t.get_where_list("name == n")
    else: # Legacy Pytables 2.x syntax
        wl = t.getWhereList("name == n")
        
    nstep = t[wl]["value"]
    if hasattr(nstep, '__len__') and 0 < len(nstep):
        nstep = nstep[0]

    n = "%-80s" % "time"
    t = h5file.getNode("/real scalars")
    if int(tables.__version__.split(".")[0]) >= 3: # Pytables 3.x syntax
        wl = t.get_where_list("name == n")
    else: # Legacy Pytables 2.x syntax
        wl = t.getWhereList("name == n")
    
    time = t[wl]["value"]
    if hasattr(time, '__len__') and 0 < len(time):
        time = time[0]

    # Check to see if the RayData dataset exists in the plot file. If
    # not, skip the file:
    if not "RayData" in h5file.root:
        print("No ray data in \"" + filename + "\" skipping...\n")
        h5file.close()
        return

    data = h5file.root.RayData[:,:]
    nrow = data.shape[0]
    h5file.close()

    # Sort the data:
    tags = data[:,0]
    indx = tags.argsort(kind="mergesort")

    sorted_data = np.empty((nrow,5))
    for i in range(len(indx)):
        sorted_data[i,:] = data[indx[i],:]

    # count = 0
    # for row in sorted_data:
    #     count += 1
    #     tag, x,y,z, power = row    
    #     print "%5i %i %15.6e %15.6e %15.6e %15.6e" % (count, tag, x,y,z, power)

    # print "\n\n"
    # print sorted_data.shape

    # Create the VTK file:
    pos = []
    powers = []
    connectivity = []

    tag = sorted_data[0,0]
    count = 0
    for i in range(nrow):
        count += 1
    
        pos.append(sorted_data[i,1])
        pos.append(sorted_data[i,2])
        pos.append(sorted_data[i,3])
        powers.append(sorted_data[i,4]* 1.0e-07) 

        if count > 1: 
            connectivity.append((3, i-1,i))

        if i == nrow-1 or tag != sorted_data[i+1,0]:
            if i != nrow-1: tag = sorted_data[i+1,0]
            count = 0

    sim = fn.split("_")[0]
    cycle = fn.split("_")[-1]
    out_fn = sim + "_las_" + str(cycle)

    # vars = (("RayPower", 1, 1, powers),)
    # visit_writer.WriteUnstructuredMesh(out_fn, 0, pos,
    #                                    connectivity, vars)

    fhand = open(out_fn + ".vtk", "w")
    fhand.write("# vtk DataFile Version 2.0\n")
    fhand.write("Written using VisIt writer\n")
    fhand.write("ASCII\n")
    fhand.write("DATASET UNSTRUCTURED_GRID\n")

    # Write out time and cycle information:
    fhand.write("FIELD FieldData 2\n")
    fhand.write("TIME 1 1 double\n")
    fhand.write(str(time) + "\n")
    fhand.write("CYCLE 1 1 int\n")
    fhand.write(str(nstep) + "\n")


    fhand.write("POINTS %i double\n" % (len(pos)/3))

    count = 0
    for i in range(len(pos)):
        fhand.write("%20.12e " % pos[i])
        count = count + 1
        if count == 9:
            fhand.write("\n")
            count = 0
    if count != 0:
        fhand.write("\n")

    fhand.write("CELLS %i %i\n" % (len(connectivity), 3*len(connectivity)))
    for i in range(len(connectivity)):
        fhand.write("2 %i %i \n" % (connectivity[i][1], connectivity[i][2]))

    fhand.write("CELL_TYPES %i\n" % len(connectivity))
    for i in range(len(connectivity)):
        fhand.write("3 \n")

    fhand.write("CELL_DATA %i\n" % len(connectivity))
    fhand.write("POINT_DATA %i\n" % (len(pos)/3))
    fhand.write("SCALARS RayPower_Watts double\n")
    fhand.write("LOOKUP_TABLE default\n")

    count = 0
    for i in range(len(powers)):
        if powers[i] < 1.0e-300: powers[i] = 0.0
        fhand.write("%20.12e " % powers[i])
        count = count + 1
        if count == 9:
            fhand.write("\n")
            count = 0
    if count != 0:
        fhand.write("\n")
    fhand.write("\n")

filenames = argv[1:]
for fn in filenames:
    write_vtk(fn)
