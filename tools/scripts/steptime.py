#!/usr/bin/env python
from optparse import OptionParser
import re
import time
from datetime import datetime
import numpy as np

def reshape(arr, cs):
    arr = np.array(arr)
    max_length = int(len(arr) / cs) * cs
    arr = arr[:max_length]
    arr.shape = (max_length/cs, cs)
    arr = arr.sum(axis=1) / cs
    return arr


def get_num_tasks(line):
    """
    Search a line to see if it contains the number of tasks. If it
    does, return it, otherwise return None
    """
    m = re.search(r"Number of MPI tasks: *([0-9]+)", line)
    if m == None: return None
    return int(m.group(1))

def get_nxb(line, direction=None):
    """
    Returns the number of cells in a block in each direction
    """

    if direction != "x" and direction != "y" and direction != "z":
        raise ValueError("Direction must be x, y, or z, got: " + repr(direction))

    m = re.search(r"Number " + direction + r" zones: *([0-9]+)", line)
    if m == None: return None
    return int(m.group(1))

def get_num_leaf_blocks(line1, line2):
    # Check if 
    if not line1.startswith(" [GRID amr_refine_derefine] min leaf blks"):
        return None
    
    line = line1
    if line2[0:2] != " [":
        line = line1 + line2

    return int(line.split()[-1])
    
def checkline(line):
    m = re.search(r".*step: n=([0-9]+).*", line)
    if m == None: return None
    return int(m.group(1))

def extracttime(line):
    m = re.match(r" \[ (\d\d-\d\d-\d\d\d\d  \d\d:\d\d:\d\d).(\d\d\d)", line)

    if m == None: return None, None, None

    date_string = m.group(1)
    msec = int(m.group(2)) * 0.001

    time_struct = datetime.strptime(date_string, "%m-%d-%Y  %H:%M:%S").timetuple()
    wall_time = time.mktime(time_struct) + msec

    m = re.search(r" t=([0-9]\.[0-9]+E[\+\-][0-9][0-9])", line)
    sim_time = float(m.group(1))

    m = re.search(r"dt=([0-9]\.[0-9]+E[\+\-][0-9][0-9])", line)
    sim_dt = float(m.group(1))

    return wall_time, sim_time, sim_dt

def main():
    usage = """usage: %prog LOGFILE

This utility reads a FLASH log file and print the amount of time in
seconds spent on each time step.

WARNING: steptime.py only works for log files containing a single run!
"""

    parser = OptionParser(usage)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
    parser.add_option("--chunk-size", action="store", type="int", dest="chunk_size")

    (opts, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error("Incorrect arguments")
    fn = args[0]
    lines = open(fn, "r").readlines()

    steps = []
    wall_times = []
    sim_times = []
    dts   = []
    dtdsim = []
    sim_dts = []
    dts_block = []
    dtdsim_block = []
    num_leaf = []

    # Load the chunk size:
    if opts.chunk_size == None:
        chunk_size = 1
    else:
        chunk_size = opts.chunk_size

    for i,line in enumerate(lines[:-1]):
        nt = get_num_tasks(line)
        if nt != None:
            num_tasks = nt

        n = get_nxb(line, "x")
        if n != None: nxb = n

        n = get_nxb(line, "y")
        if n != None: nyb = n

        n = get_nxb(line, "z")
        if n != None: nzb = n

        # Get number of leaf blocks:
        nlb_new = get_num_leaf_blocks(lines[i], lines[i+1])

        if nlb_new != None:
            nlb = nlb_new
            
        step = checkline(line)
        if step == None: continue

        # This is a line containing time step information:
        wall_time, sim_time, sim_dt = extracttime(line)
        if wall_time == None: 
            raise ValueError("step %d line %d has no time stamp" % (step,i))

        steps.append(step)
        wall_times.append(wall_time)
        sim_times.append(sim_time)
        sim_dts.append(sim_dt)
        num_leaf.append(nlb)

    for i in xrange(len(steps)-1):
        dts.append(wall_times[i+1]-wall_times[i])
        dtdsim.append(dts[-1]/sim_dts[i])
        dts_block.append(dts[-1]/num_leaf[i])
        dtdsim_block.append(dtdsim[-1]/num_leaf[i])
        

    # Compute chunked time step length:
    dts = reshape(dts, chunk_size)
    dtdsim = reshape(dtdsim, chunk_size)
    num_leaf = reshape(num_leaf, chunk_size)
    sim_times = reshape(sim_times, chunk_size)
    sim_dts = reshape(sim_dts, chunk_size)
    

    print "# Num Tasks  =", num_tasks
    print "# Chunk Size =", chunk_size
    print "# NXB =", nxb
    print "# NYB =", nyb
    print "# NZB =", nzb
    print ""
    print "# Columns are:"
    print "%6s %20s %20s %10s %20s %20s %20s %20s %20s %20s" % \
        ( "# step", 
          "sim time",
          "dt",
          "leaf blks", 
          "wall",
          "wall/dt",
          "wall/cell",
          "wall/dt/cell",
          "CPU-hrs/cell/step",
          "wall*proc/dt/cell" )
    for i in xrange(len(dts)):
        print "%6d %20.6e %20.6e %10i %20.6e %20.6e %20.6e %20.6e %20.6e %20.6e" % \
            ( steps[i*chunk_size],
              sim_times[i],
              sim_dts[i],
              int(num_leaf[i]),
              dts[i],
              dtdsim[i],
              dts_block[i]/(nxb*nyb*nzb),
              dtdsim_block[i]/(nxb*nyb*nzb),
              dts_block[i]*num_tasks / (nxb*nyb*nzb) / 3600.0,
              dtdsim_block[i] * num_tasks / (nxb*nyb*nzb))

if __name__ == "__main__":
    main()


