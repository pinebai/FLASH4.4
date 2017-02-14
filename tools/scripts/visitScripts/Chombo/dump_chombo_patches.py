import os, sys, glob

chkdir = os.getcwd()
chkglob = "*_hdf5_chk_[0-9][0-9][0-9][0-9]"
fullglob = chkdir + "/" + chkglob
chknames = glob.glob(fullglob)
chknames.sort()
print "The checkpoint names are:", chknames

#Here are all the specific visit commands.
for chkname in chknames:
  SaveWindowAtts = SaveWindowAttributes()
  SaveWindowAtts.format = SaveWindowAtts.PNG
  SaveWindowAtts.fileName = chkname
  SaveWindowAtts.family = 0
  SetSaveWindowAttributes(SaveWindowAtts)

  OpenDatabase(chkname, 0)
  AddPlot("Mesh", "Mesh", 1, 1)
  AddPlot("Pseudocolor", "dens", 1, 1)
  AddPlot("Subset", "patches", 1, 1)
  ChangeActivePlotsVar("patches")
  SubsetAtts = SubsetAttributes()
  SubsetAtts.singleColor = (0, 0, 0, 255)
  SubsetAtts.colorType = SubsetAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
  SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
  SubsetAtts.lineWidth = 4
  SubsetAtts.wireframe = 1
  SetPlotOptions(SubsetAtts)
  DrawPlots()
  SaveWindow()
  DeleteAllPlots()
  CloseDatabase(chkname)

CloseComputeEngine()
sys.exit()
