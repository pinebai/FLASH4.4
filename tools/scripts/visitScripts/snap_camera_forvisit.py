#!/usr/bin/env python

#================================================================================
# VARIABLES YOU MAY CHANGE (see below "import" statements:
# camvec:  vector to camera (presumably from somewhere near center of data)
#    default:  (0, 0, 1)
#    (referred to in "VisIt User's Manual" as "view normal", but shown as a vector
#    FROM the camera TO the camera focus)
# zoom:  factor by which view is zoomed into the image
#      1.0:  not zoomed
#    < 1.0:  zoomed out
#    > 1.0:  zoomed in
# reflev:  refinement level to visualize (from 0 to n-1)
#    can be used for fast, low-res rendering
#    There's no programmatic way to determine the valid level numbers;
#    the user needs to know them
# viewsscope_flag:  choose whether view includes entire data domain or just the plot
#    0:  entire data domain
#    1:  only the plot
# varname:  name of scalar variable to visualize
# plotfilename:  name of input file
# imagefilename:  base of name PNG output image
#================================================================================
#================================================================================
# OTHER VARIABLES:
# upvec:  up vector for graphical viewing process
#    default:  (0, 1, 0)
#    changed by script if camvec is also (0, 1, 0)
#================================================================================

import sys
import os

#================================================================================
# YOU MAY CHANGE THESE
#================================================================================
varname = "dens"
#plotfilename = "sedov_hdf5_plt_cnt_0008"
plotfilename = "c_hdf5_plt_cnt_0105"
imagefilename = "RT_dens_isosurface_onlyplot"
camvec = (0.5, 0.5, 1.0)
zoom = 0.75
reflev = 1
viewscope_flag = 1
#================================================================================

upvec = []
upvec_plus_y = [0.0, 1.0, 0.0]
upvec_minus_z = [0.0, 0.0, -1.0]

# Plot-snapper
def snap():
	 global plotfilename, varname, imagefilename
	 global upvec, camvec, zoom, reflev

	 OpenDatabase(plotfilename)
	 if viewscope_flag == 0:
		  SetViewExtentsType("original") # Show entire data dmoain
	 elif viewscope_flag == 1:
		  SetViewExtentsType("actual") # Crops view tightly about the plot
	 availableplots = PlotPlugins()
#	 print "availableplots: ", availableplots

#	 if ("Pseudocolor" in availableplots):
#		  AddPlot("Pseudocolor", varname)
	 if ("Contour" in availableplots):
		  AddPlot("Contour", varname)
		  conatts = ContourAttributes()
		  conatts.contourNLevels = 1 # reduce num contours to a zippier '1'
		  SetPlotOptions(conatts)

		  AddOperator("IndexSelect") # "Index select" op. to select lo-res level for fast rendering
		  idxsel_atts = IndexSelectAttributes()
		  idxsel_atts.groupIndex = reflev
		  SetOperatorOptions(idxsel_atts)

		  DrawPlots()

		  view3d = GetView3D()
#		  print "The given 3d view is: ", view3d

		  view3d.viewNormal = camvec
		  view3d.viewUp = (tuple(upvec))
		  view3d.imageZoom = zoom
#		  print "The changed 3d view is: ", view3d
		  SetView3D(view3d)

		  swa = SaveWindowAttributes()
		  swa.format = swa.PNG
		  swa.fileName = imagefilename
		  swa.family = 0 # not saving series ("family") of images
		  SetSaveWindowAttributes(swa)
		  SaveWindow()
	 else:
		  print "Contour plot not available"
		  print "Exiting..."
		  sys.exit()

#================================================================================
# RUN
#================================================================================

if camvec == tuple(upvec_plus_y):
	 upvec = tuple(upvec_minus_z)
	 print "[main] - camvec == upvec_plus_y, so upvec will NOT be default: ", upvec
else:
	 upvec = tuple(upvec_plus_y)
	 print "[main] - camvec != upvec_plus_y, so upvec WILL be default: ", upvec

snap()
sys.exit()

