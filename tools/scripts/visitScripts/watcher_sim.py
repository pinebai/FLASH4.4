#!/usr/bin/env python

# This script simulates the calling, by Bob Fisher's "watcher" script,
# of snap3d_viewpoint_pseudocolor.py

import os
import glob
import sys

filelist = []
#filelist = glob.glob('c*plt*')
#filelist = glob.glob('sedov_hdf5_plt_cnt_0008')
#filelist = glob.glob('gw_3D6_7_hdf5_plt_crn_0050')
filelist = glob.glob('c_hdf5_plt_cnt_0105')

scriptname = "/Users/hudson/Projects/Flash3/tools/scripts_WorkingCopy/snap_camera_forwatcher.py"
varname = "dens"
camvec = "1 1 1"

i = 0
for pfname in filelist:
	 iname = "imageout" + "%04d"%(i,)
	 os.system("visit -cli -nowin -default_format FLASH -s %s --pfname %s --vname %s --ifname %s --camvec %s" % (scriptname, pfname, varname, iname, camvec))
#	 os.system("visit -cli -default_format FLASH -s %s --pfname %s --vname %s --ifname %s --camvec %s" % (scriptname, pfname, varname, iname, camvec))
	 i += 1

