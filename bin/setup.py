#!/usr/bin/env python

import os, sys, string, re, time, shutil, UserDict, types, glob, socket, math

import parseCmd, lazyFile
import globals
from globals import *   # GVars and SetupError
from utils import *     # Assorted little cute functions
from libUtils import *  # for libUnion
from linkFiles import * # for LinkFileList class
from genFiles import *  # code to generate files (Makefiles, Flash.h, etc)
from unitUtils import * # for UnitList class

########################### START SCRIPT ################################
def main():
    sys.setcheckinterval(10000) # this is not threaded application
 
    parseCmd.init()
    # setup.py is FLASH_HOME/bin/setup.py
    # FLASH_HOME is the parent of directory having setup.py 
    flashHome = os.path.dirname(
                   os.path.dirname(os.path.abspath(sys.argv[0])))
    # Initialize other related directory names
    GVars.init(flashHome)

    # parse the command line arguments. The effect of calling
    # this method is to set a lot of values in 'GVars'
    parseCmd.parse()

    ############## Clean up simulation name ############################

    os.chdir(GVars.simulationsDir)
    if not os.path.isdir(GVars.simulationName): # Needs to be cleaned up?
       # Find directories matching simulatation name (case insensitive)
       simnames = dirGlob(GVars.simulationName+"*") 

       if len(simnames) == 1: # only one possible name
          GVars.simulationName = simnames[0]
       else:
          if not simnames: # no match
             raise SetupError("Unit %s not found" % GVars.simulationName)
          else: 
             GVars.out.put("Potential completions : %s"%', '.join(simnames),globals.ERROR)
             raise SetupError("Unit %s is ambiguous" % GVars.simulationName)

    ############ Proceed with other stuff ###################	 
    os.chdir(GVars.flashHomeDir)

    machDir = determineMachine()

    unitList = UnitList() # A class which encapsulates operations of unit collections

    objdir = os.path.join(GVars.flashHomeDir,GVars.objectDir)
    if not os.path.isdir(objdir): os.mkdir(objdir)

    if GVars.unitsFile: # copy over the user specified units file
       sname = os.path.join(GVars.flashHomeDir,GVars.unitsFile)
       tname = os.path.join(objdir,globals.UnitsFilename)
       shutil.copy(sname,tname)

    ############## SWITCHING TO SOURCE DIR ############################
    os.chdir(GVars.sourceDir)

    if GVars.auto:
        # generate and write an automatic "Units" file
	unitList.generateUnitsfile()

    # read units file, satisfy requirement and all the good stuf
    unitList.populate() 

    # adjust values of command line options based on units used
    # this is the only place the setup script uses its knowledge
    # about Flash Units and what they do
    unitList.adjustOpts()
    parseCmd.final() # finalise the options

    # create the object which does stuff relating to linking files
    linkList = LinkFileList(objdir)

    # Combine info from all Config files 
    configInfo = unitList.getConfigInfo() 
    # get Runtime Paramas info
    rpInfo = unitList.getRPInfo(configInfo['max_plot_vars'])
    # get info regarding all variables
    varInfo = unitList.getVarInfo()

    # Add library options to configInfo
    for key in GVars.withLibraries.keys():
        if not configInfo['LIBRARY'].has_key(key): configInfo['LIBRARY'][key] = []
        configInfo['LIBRARY'][key].append(GVars.withLibraries[key])

    ######### Now handle libraries depending on other libraries ##########
    configInfo['libConfigInfo'] = LibUnion(configInfo['LIBRARY'])
    # ConfigInfo['libConfigInfo'] is a dictionary mapping libraries to their config information
    # also CLASS.libOrder is a list of pairs (libname, args)
    # This order is important due to linker invocation order 

    configInfo['libConfigInfo'].writeLibraries() # write out lib data to setup_libraries

    # only now is noClobber flag stable
    linkList.cleanObjectDir()
    ############################### CHANGE TO OBJECT DIR ###################
    os.chdir(GVars.flashHomeDir)
    os.chdir(GVars.objectDir)
    # remove the SuccessFile
    try:
       os.unlink(globals.SuccessFilename)
    except OSError: 
       pass

    rpInfo.writeCode(configInfo) # write Fortran code for RP support

    
    # find files which should not be linked
    linkList.getDontLinkList(unitList.getList(),configInfo['LINKIF'])

    # getLinkOrder does the fancy sorting of unitnames
    for unitname in unitList.getLinkOrder():
        linkList.linkFiles(os.path.join(GVars.sourceDir, unitname))
    linkList.linkFiles(os.path.join(GVars.flashHomeDir, 'sites'))

    # Link in the right version of the make file
    # attempt to find Makefile in FLASH root directory, otherwise link make file from sites
    tmp_makefilePathRoot=os.path.join(flashHome, "Makefile.h" + GVars.makefileext)
    tmp_makefilePathSites= os.path.join(machDir, "Makefile.h"+GVars.makefileext)

    if os.path.isfile(tmp_makefilePathRoot):
    	linkList.addLink(tmp_makefilePathRoot, "Makefile.h")
        print "Using Makefile.h: " + tmp_makefilePathRoot	
    else:
	linkList.addLink(tmp_makefilePathSites,"Makefile.h") 
	print "Using Makefile.h: " + tmp_makefilePathSites

    # functions in the simulations dir override earlier instances
    linkList.linkFiles(GVars.simulationsDir)

    # functions from LINKIF statements override everything else,
    # assuming their conditions apply.
    linkList.doLINKIFOverrides(unitList.getList(), configInfo['LINKIF'])

    # now is when we do the real link/copying
    linkList.reallyLink()

    ############## flash.par and Makefiles **************

    #  Copy in flash.par
    if os.path.isfile(os.path.join(GVars.simulationsDir, GVars.simulationName, 'flash.par')):
        shutil.copy(os.path.join(GVars.simulationsDir, GVars.simulationName, 'flash.par'), '.')

    test_src = os.path.join(GVars.flashHomeDir, 'tools', 'scripts', 'testing',
                            'precision_test', 'precision_test.F90')    
    if os.path.isfile(test_src):
        shutil.copy(test_src, '.')

    GVars.out.put('creating Makefiles for all units',globals.INFO)

    unitList.createMakefiles()
    #FIXME merge with createMakefiles

    GVars.out.put('generating buildstamp generator',globals.INFO)
    generateBuildstampGenerator()

    GVars.out.put('copying release accessor function Makefile',globals.INFO)
    shutil.copy(os.path.join(GVars.flashHomeDir, 'bin/make_release'), '.')

    GVars.out.put('copying buildstats accessor function Makefile',globals.INFO)
    shutil.copy(os.path.join(GVars.flashHomeDir, 'bin/make_bstats'), '.')

    GVars.out.put('copying flashUnits accessor function Makefile',globals.INFO)
    unitList.generateSetupFlashUnits()
    #shutil.copy(os.path.join(GVars.flashHomeDir, 'bin/make_flashUnits'), '.')

    shutil.copy(os.path.join(GVars.flashHomeDir, 'bin/resetup'), '.')

    GVars.out.put('copying Dependency generator',globals.INFO)
    flist = ["setup_depends.py","setup_addcdepends.py","setup_reorder.py","reorder.tpl"]
    for f in flist:
        shutil.copy(os.path.join(GVars.flashHomeDir, 'bin/'+f), '.')

    GVars.out.put('generating Makefile',globals.IMPINFO)

    generateMakefile(configInfo, machDir)
    
#    generateDbaseDefines(configInfo, opts)

    generateFlashDefines(configInfo)  # writes the Flash.h file
    writeSimulationFiles(configInfo)
   
    # Copy/link datafiles to object directory
    if GVars.portable:
      func = lambda real,link : shutil.copy2(real,link)
    else:
      cwd = os.getcwd()
      func = lambda real,link : os.symlink(getRelPath(real,cwd),link)

    datafiles = []
    for wildcard in configInfo['DATAFILES']:
        for file in glob.glob(os.path.join(GVars.sourceDir,wildcard)):
            datafiles.append(file[len(GVars.sourceDir):].strip("/"))
            bname = os.path.basename(file)
            try:
               func(file,bname)
            except:
               os.remove(bname)
               func(file,bname)
    simDir = os.path.join(GVars.simulationsDir, GVars.simulationName)
    for wildcard in GVars.datafiles:
        for file in glob.glob(os.path.join(simDir,wildcard)):
            datafiles.append(file[len(GVars.sourceDir):].strip("/"))
            bname = os.path.basename(file)
            try:
               func(file,bname)
            except:
               os.remove(bname)
               func(file,bname)
    if datafiles: 
       GVars.out.put('Copying data files: %d copied' % len(datafiles),globals.IMPINFO)
    # write "setup_datafiles" into object directory, even if empty
    open(os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.SetupDatafilesFilename), 'w').write("\n".join(datafiles))
    if (GVars.setupVars.get("ParameshLibraryMode") or 
        ((not GVars.setupVars.get("Grid") or (GVars.setupVars.get("Grid").upper()=="PM4DEV")) and
         (configInfo['PPDEFINES'].has_key('USE_AMR_RUNTIME_PARAMETERS_FILE') or
          ('-DUSE_AMR_RUNTIME_PARAMETERS_FILE' in GVars.definesNames)) )):
        generatePm3RPDefines(configInfo)
    if GVars.parfile:
       shutil.copy2(os.path.join(simDir,GVars.parfile),"flash.par")
       GVars.out.put("Copied %s as flash.par" % GVars.parfile,globals.IMPINFO)

    parseCmd.writeCmdLine()
    unitList.writeSetupUnitsFile()
    varInfo.writeVarInfo()
    # rp Stuff
    rpInfo.writeRPInfo()
    rpInfo.writeDefaultPar()
    # if no flash.par copy default.par over
    if not os.path.isfile('flash.par'):
       shutil.copy(globals.RPDefaultParFilename, 'flash.par')
   
    # create the successfile
    ofd = open(globals.SuccessFilename,"w")
    ofd.write("SUCCESS\n")
    ofd.close()

    os.chdir(GVars.flashHomeDir)
    GVars.out.put('SUCCESS',globals.ALWAYS)

########################### END SCRIPT #######################################
if __name__=='__main__':
    try:
        main()
    except SetupError, inst:
        if inst.args: print inst.args[0]
        sys.exit(1)
    except KeyboardInterrupt:
        print '\nuser abort'
        sys.exit(2)
    except Exception, e:
        if hasattr(e, 'dumbuser'):
            print '\nUser Error!!!\n%s' % e.message
        else:
            print '\nA setup internal error has occured, if possible please email the following\ndebugging info to flash-bugs@flash.uchicago.edu'
            print 'Arguments:', sys.argv
            print 'Python Version: %d.%d.%d' % sys.version_info[:3]
            print 'Platform Details: %s' % sys.platform
            raise e

