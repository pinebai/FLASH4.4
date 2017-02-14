# #!/usr/bin/env python

#================================================================================
# VARIABLES YOU MAY CHANGE:
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
import optparse

zoom = 0.75

reflev = 1

viewscope_flag = 1

upvec = []
upvec_plus_y = [0.0, 1.0, 0.0]
upvec_minus_z = [0.0, 0.0, -1.0]

#--------------------------------------------------------------------------------
# Plot-snapper
def snap():
     global opt
     global upvec, zoom, reflev

     OpenDatabase(opt.plotfilename)
     if viewscope_flag == 0:
          SetViewExtentsType("original") # Show entire data domain
     elif viewscope_flag == 1:
          SetViewExtentsType("actual") # Crops view tightly about the plot
     availableplots = PlotPlugins()
     print "availableplots: ", availableplots

#	 if GetGlobalAttributes().autoUpdateFlag != 0:
#		  print "Auto update flag is ", GetGlobalAttributes().autoUpdateFlag
#		  print "So, the necessary, initial creation of a high-res plot is (expensively) rendered"

#	 if ("Pseudocolor" in availableplots):
#		  AddPlot("Pseudocolor", opt.varname)
     if ("Contour" in availableplots):
          AddPlot("Contour", opt.varname)
#		  DrawPlots()
          con_atts = ContourAttributes()
          con_atts.contourNLevels = 1
          SetPlotOptions(con_atts)

          AddOperator("IndexSelect") # "Index select" op. to select lo-res level for fast rendering
          idxsel_atts = IndexSelectAttributes()
          idxsel_atts.groupIndex = reflev
          SetOperatorOptions(idxsel_atts)

          DrawPlots()

          view3d = GetView3D()
          print "The given 3d view is: ", view3d

          view3d.viewNormal = opt.camvec
          view3d.viewUp = (tuple(upvec))
          view3d.imageZoom = zoom
		  print "The changed 3d view is: ", view3d
		  SetView3D(view3d)

		  swa = SaveWindowAttributes()
		  swa.format = swa.PNG
		  swa.fileName = opt.imagefilename
		  swa.family = 0 # not saving series ("family") of images
		  SetSaveWindowAttributes(swa)
		  SaveWindow()
	 else:
		  print "Contour plot not available"
		  print "Exiting..."
		  sys.exit()

#================================================================================
# RUN
# Let's try the "optparse" module
#    (http://www.python.org/doc/2.3.5/lib/module-optparse.html)
#    (http://www.python.org/doc/2.3.5/lib/optparse-philosophy.html)
#================================================================================
print "Length of sys.argv is %s" % len(sys.argv)

parser = optparse.OptionParser()

# print "[main] - argv: ", sys.argv

# To get visit args off my back (to get it past optparse's persnickety restrictions:
#parser.add_option("-cli")
parser.add_option("-d") # for "-default_format"
parser.add_option("-s") # for "-s"
parser.add_option("-n") # for "-nowin"

# For this script's args:
parser.add_option("--pfname", action="store", type="string", dest="plotfilename")
parser.add_option("--ifname", action="store", type="string", dest="imagefilename")
parser.add_option("--vname", action="store", type="string", dest="varname")
parser.add_option("--camvec", action="store", type="float", dest="camvec", nargs=3)

parser.set_defaults(camvec=(0, 0, 1))

(opt, args) = parser.parse_args()
# print "[main] - args                   : ", args
# print "[main] - opt                    : ", opt

if opt.imagefilename == None:
	 print "---------- imagefilename is None ------------"
	 print "Exiting..."
	 sys.exit()

if opt.plotfilename == None:
	 print "---------- plotfilename is None ------------"
	 print "Exiting..."
	 sys.exit()

if opt.varname == None:
	 print "---------- varname is None ------------"
	 print "Exiting..."
	 sys.exit()

print "[main] - plot file name: ", opt.plotfilename
print "[main] - image file name: ", opt.imagefilename
print "[main] - variable name: ", opt.varname
print "[main] - direction-of-camera vector: ", opt.camvec
print "[main] - args                   : ", args

if opt.camvec == tuple(upvec_plus_y):
	 upvec = tuple(upvec_minus_z)
	 print "[main] - camvec == upvec_plus_y, so upvec will NOT be default: ", upvec
else:
	 upvec = tuple(upvec_plus_y)
	 print "[main] - camvec != upvec_plus_y, so upvec WILL be default: ", upvec

# print "[main] - up vector as list: ", upvec

snap()
sys.exit()

