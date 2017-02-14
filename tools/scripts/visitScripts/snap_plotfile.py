

#================================================================================
# README
#================================================================================


#================================================================================
# DEFINITIONS
#================================================================================
import sys
import os

#varname = []
#plotfilename = []
#imagefilename = []
varname = ""
plotfilename = ""
imagefilename = ""

# Handlers of command-line arguments
def handle_plotfilename(pname):
	 global plotfilename
	 plotfilename = pname
	 print "Plot file name is '%s'" % (plotfilename,)
def handle_variablename(vname):
	 global varname
	 varname = vname
	 print "Variable name is '%s'" % (varname,)
def handle_imagename(iname):
	 global imagefilename
#	 imagefilename = iname + ".png"
	 imagefilename = iname
	 print "Image file name is '%s'" % (imagefilename,)

# Dictionary of functions for valid command-line keywords
functionDict = {
	 "-fn": handle_plotfilename,
	 "-vn": handle_variablename,
	 "-in": handle_imagename,
	 }

# Dictionary for command-line arguments
argumentDict = {
	 }

# Plot-snapper
def snap():
	 global plotfilename, varname, imagefilename
	 OpenDatabase(plotfilename)
	 availableplots = PlotPlugins()
	 availableops = OperatorPlugins()
#	 for anop in availableops:
#		  print anop
#	 sys.exit()
	 if ("Pseudocolor" in availableplots) and ("Slice" in availableops):
		  AddPlot("Pseudocolor", varname)
		  AddOperator("Slice")
		  sa = SliceAttributes()
		  sa.project2d = 0
		  sa.normal = 0.0, 0.0, 1.0
		  sa.originType = 2
		  sa.originPercent = 0
		  SetOperatorOptions(sa)
		  DrawPlots()
		  swa = SaveWindowAttributes()
		  swa.format = swa.PNG
		  swa.fileName = imagefilename
		  swa.family = 0 # not saving series ("family") of images
		  SetSaveWindowAttributes(swa)
		  SaveWindow()
	 else:
		  print "At least one of {Pseudocolor plot, Slice operator} not available"
		  print "Exiting..."
		  sys.exit()

#================================================================================
# RUN
#================================================================================
print "Length of sys.argv is %s" % len(sys.argv)

if len(sys.argv) == 11:
	 for i in range(5, 10, 2):
		  print "argv[%d]: %s" % (i, sys.argv[i])
		  print "argv[%d]: %s" % (i+1, sys.argv[i+1])
		  argumentDict[sys.argv[i]] = sys.argv[i+1]
		  print "argumentDict[%s]: %s" % (sys.argv[i], argumentDict[sys.argv[i]])
	 print "Dictionary of command-line arg's built."
	 for keyword in functionDict.keys():
		  print "Keyword from functionDict: %s" % (keyword,)
		  if argumentDict.has_key(keyword):
				functionDict[keyword](argumentDict[keyword])
		  else:
				print "%s is an invalid command-line keyword" % (keyword,)
				print "Exiting script..."
				sys.exit()
	 snap()
	 sys.exit()
else:
	 print "There must be three argument pairs:  '-fn <plotfilename>', '-vn <varname>' and '-in <imagefilename>'"
	 print "(The pairs can appear in any order)"
	 print "This script did not create a plot or an image thereof"
	 sys.exit()

