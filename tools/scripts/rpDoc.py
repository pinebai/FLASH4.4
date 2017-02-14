#!/usr/bin/env python2.7
import os.path, sys, re

# rpDoc.py is FLASH_HOME/tools/scripts/setup.py
# FLASH_HOME is the parent of directory having setup.py 
flashHome = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))
sys.path.insert(0, os.path.join(flashHome, "bin"))

from globals import *   # GVars, opts and SetupError
from unitUtils import * # for UnitList class

unitreg = re.compile(r"(^|/)[A-Z][^/]*(?=(/|$))")

def pullunitname(name):
    m = unitreg.search(name)
    if not m:
       print name," did not match. Programming Error"
       sys.exit(1)
    else:
       return name[:m.end()]

########################### START SCRIPT ################################
def main():
   # Initialize other related directory names
   GVars.init(flashHome)
   os.chdir(GVars.sourceDir)

   ulist = [] # list of all units
   for root, dirs, files in os.walk("."):
     if ".svn" in dirs:
       dirs.remove(".svn")
     dirs.sort()
     ulist.extend([os.path.normpath(os.path.join(root, dir)) for dir in dirs])

   unitList = UnitList(ignorePP=True)
   unitList.addUnits(ulist)
   # The following string argument will be appended to "plot_var_" to
   # create an entry for "plot_var_<N>, ..." in the report files.
   # Replace by 0 if you prefer.
   rpInfo = unitList.getRPInfo("<N>, for N=1..MAX_PLOT_VARS",warn=None)

   # list of top level locations
   locs = set( [ pullunitname(x[0]) for x in rpInfo.keys() ] )
   ddocs = os.path.join(GVars.flashHomeDir,"docs","designDocs")
   rpInfo.writeRPInfo(os.path.join(ddocs,"rpDoc.txt"))
   rpInfo.writeDuplications(os.path.join(ddocs,"rpDuplications.txt"))
   for ul in locs:
       lastpart = ul.split("/")[-1]
       rpInfo.writeRPInfo(os.path.join(ddocs,"rp_%s.txt"%lastpart),prefix=ul)


########################### END SCRIPT #######################################
if __name__=='__main__':
    main()

