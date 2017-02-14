#!/usr/bin/env python
import os.path, sys, string, re, getopt
try:
   import subprocess  # use this when available to avoid a DeprecationWarning
except ImportError:
   import popen2
from interfacesDict import *  # for InterfacesDict class and associated constants

########################### codeCheck class ############################

class checkCode(dict):
   """This class checks a file for coding violations. To add a new 
check add the name of the new violation to the "violations" variable, 
and create methods <violation>_<ext> for each extension for which this
violation applies. This method takes as input a list of strings, 
representing the contents of the file and should return some 
datastructure representing the violations. Usually it is just the list 
of line numbers where the violation occurred. The structure of this 
datastructure should depend only the violation in question.

If the violation does not apply for a certain extension,
you dont have to create the method (or create one which always 
returns [])

Also for each violation you need to create a method <violation>_report 
which takes a list of pairs of (filename, appropriate datastructure) 
and the filename which should contain the report and generate the report 
containing the violations. 

The class also implements fixes. To add a new fix, add a dictionary 
entry to the fixes list, and create a method called fix_<nickname>_<ext> 
where <nickname> is the nickname you gave for the fix, and ext is the 
fileextension to which it should apply. If you define one without the 
_<ext> it will apply to files with all extensions. The method to fix, 
will get as input a list of strings representing the file contents.

It should modify the list in place. The return value should be
None -- if no change was made to the list 1 -- if some change was made 
to the list (so the file will be updated on disk)

The violation/fix dictionaries also support these optional entries:

"exceptions" : if present, should be keyed to a list of regular
expressions by which directories may be excluded from a check. The
violation or fix will not be applied if the directory/file matches a
pattern in the "exceptions" list.

"isComplex" : if present, should be keyed to a method which accepts
a single argument "lines", which represents a list of the lines of
the file being read. If, upon examining "lines" the method returns
True, the file is considered "too complex" and will not be checked
for the violation. The file will appear in a list of complex files
at the beginning of the report.
If no "isComplex" method is provided, this script will use the default
hard-coded "isComplex" method.

"disable": if present, should be keyed to "True". Setting "disable"
simply prevents this violation / fix from being applied to any file.
"""

   pmTwo = "source/Grid/GridMain/paramesh/Paramesh2"
   pmFour0Kernel = "source/Grid/GridMain/paramesh/paramesh4/Paramesh4.0/PM4_package"
   pmFourAny = "source/Grid/GridMain/paramesh/paramesh4/"
   pmFourDev = "source/Grid/GridMain/paramesh/paramesh4/Paramesh4dev/"
   pmFourDevKernel = "source/Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package"
   pmInterp = "source/Grid/GridMain/paramesh/interpolation/"
   PPM   = "source/physics/Hydro/HydroMain/split/PPM/PPMKernel"
   PPMexp = "source/physics/Hydro/HydroMain/split/PPM_exp"
   SPPM   = "source/physics/Hydro/HydroMain/split/SPPM/"
   SlowFlow   = "source/physics/Hydro/HydroMain/split/SlowFlowPPM/"
   IncompNS = "source/physics/IncompNS"
   IncompNSVarDens = "source/physics/IncompNS/IncompNSMain/vardens"
   Samrai = "source/Grid/GridMain/Samrai"
   MPIHybrid = "source/IO/IOMain/MPIHybrid"
   SamraiIO = "source/IO/IOMain/hdf5/parallel/Samrai"
   MHD = "source/physics/MHD/MHDMain/"
   ST = "source/physics/sourceTerms/"
   ST_B = "source/physics/sourceTerms/Burn/"
   ST_C = "source/physics/sourceTerms/Cool/"
   ST_Delep = "source/physics/sourceTerms/Deleptonize/DeleptonizeMain"
   ST_DelepLAPI = "source/physics/sourceTerms/Deleptonize/localAPI"
   ST_H = "source/physics/sourceTerms/Heat/HeatMain"
   ST_I = "source/physics/sourceTerms/Ioniz/"
   ST_ED = "source/physics/sourceTerms/EnergyDeposition/"
   DiffUG = "source/physics/Diffuse/DiffuseMain/UG"
   DiffCG = "source/physics/Diffuse/DiffuseMain/CG"
   ImBound = "source/physics/ImBound/ImBoundMain"
   GP = "source/Grid/GridParticles/"
   Part = "source/Particles"
   PartHyb = "source/Particles/ParticlesMain/active/charged/HybridPIC/"
   PartPredCorOld = "source/Particles/ParticlesMain/passive/PredictorCorrectorOld"
   IOPnet = "source/IO/IOMain/pnetcdf"
   IOdir   = "source/IO/IOMain/direct"
   IOPart = "source/IO/IOParticles"	
   Grav = "source/physics/Gravity"
   MultigridFFT = "source/Grid/GridSolvers/Multigrid/fft"
   MultigridExperimental = "source/Grid/GridSolvers/Multigrid_experimental"
   MultigridMC = "source/Grid/GridSolvers/MultigridMC"
   MultigridMC2 = "source/Grid/GridSolvers/MultigridMC_VarDens_HYPRE"
   MultigridDiffHgADI = "source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI"
   MultigridDiffHgFFT = "source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgFFT"
   paraBurn = "source/physics/sourceTerms/Burn/BurnMain/parametricBurn"
   paraBurnOld = "source/physics/sourceTerms/Burn/BurnMain/parametricBurn/Co2nse"
   linalg = "source/physics/utilities/solvers/LinearAlgebra"
   cosmoKernel = "source/physics/Cosmology/CosmologyMain/MatterLambdaKernel"
   materialProperties = "source/physics/materialProperties"
   Interfaces = "source/physics/Interfaces"
   hgFlash2 = "source/Grid/GridSolvers/Multigrid/hg_flash2"
   hgFlash2fft = "source/Grid/GridSolvers/Multigrid/hg_flash2/fft"
   BarnesHutTree = "source/Grid/GridSolvers/BHTree"
   BIPCGStab = "source/Grid/GridSolvers/BiPCGStab"
   BIPCGStabVarcoeff = "source/Grid/GridSolvers/BiPCGStab/varcoeffpoisson"
   GridStructures = "source/Grid/GridStructures"
   Aprox13t = "source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13t"
   Aprox13tInt = "source/physics/sourceTerms/Burn/BurnIntegrate/Aprox13t"
   DrIncompNS = "source/Driver/DriverMain/INSfracstep"
   EosNuc     = "source/physics/Eos/EosNuclear/kernel"
   SolidMechanics = "source/physics/SolidMechanics/SolidMechanicsMain"
   SimAdvect = "source/Simulation/SimulationMain/Advect"
   SimAdvectMass = "source/Simulation/SimulationMain/AdvectMassScalars"
   SimBlast2 = "source/Simulation/SimulationMain/Blast2"
   SimBubLab = "source/Simulation/SimulationMain/BubLab"
   SimCDMHalo = "source/Simulation/SimulationMain/CDMHalo"
   SimCell = "source/Simulation/SimulationMain/Cellular"
   SimCellpb = "source/Simulation/SimulationMain/CellularParametricBurn"
   SimCollShock = "source/Simulation/SimulationMain/CollShock"
   SimCoreColl  = "source/Simulation/SimulationMain/CoreCollapse"
   SimConduction = "source/Simulation/SimulationMain/ConductionDelta"
   SimDust = "source/Simulation/SimulationMain/DustCollapse"
   SimDustCart = "source/Simulation/SimulationMain/DustCollapseCart"
   SimFlame1 = "source/Simulation/SimulationMain/Flame1StageNoise"
   SimFlame3 = "source/Simulation/SimulationMain/Flame3StageNoise"
   SimFlameBub = "source/Simulation/SimulationMain/FlameBubble"
   SimFlash2Convert = "source/Simulation/SimulationMain/Flash2Convert"
   SimGadget = "source/Simulation/SimulationMain/GadgetSnapshot"
   SimGalaxy = "source/Simulation/SimulationMain/GalaxyCluster"
   SimGalaxyM = "source/Simulation/SimulationMain/magnetoHD/GalaxyClusterMagnetic"
   SimClustSlos   = "source/Simulation/SimulationMain/ClusterSloshing"
   SimClustSlosM  = "source/Simulation/SimulationMain/magnetoHD/ClusterSloshingMagnetic"
   SimClustSlosMN = "source/Simulation/SimulationMain/magnetoHD/ClusterSloshingMagneticNew"
   SimIV = "source/Simulation/SimulationMain/IsentropicVortex"
   SimINS = "source/Simulation/SimulationMain/INavierStokes"
   SimJeans = "source/Simulation/SimulationMain/Jeans"
   SimLayer3 = "source/Simulation/SimulationMain/Layer3"
   SimNucOToRT = "source/Simulation/SimulationMain/NucOToRT"
   SimNuc2Grid = "source/Simulation/SimulationMain/Nuc2Grid"
   SimOrbit = "source/Simulation/SimulationMain/Orbit"
   SimPancake = "source/Simulation/SimulationMain/Pancake"
   SimMacLarin = "source/Simulation/SimulationMain/MacLaurin"
   SimMHDAlfven = "source/Simulation/SimulationMain/magnetoHD/AlfvenWave"
   SimMHDBlastBS = "source/Simulation/SimulationMain/magnetoHD/BlastBS"
   SimMHDBlastMHD = "source/Simulation/SimulationMain/magnetoHD/BlastMHD"
   SimMHDBrio = "source/Simulation/SimulationMain/magnetoHD/BrioWu"
   SimMHDCloud = "source/Simulation/SimulationMain/magnetoHD/CloudShock"
   SimMHDCurrent = "source/Simulation/SimulationMain/magnetoHD/CurrentSheet"
   SimMHDDegen = "source/Simulation/SimulationMain/magnetoHD/DegenEOS"
   SimMHDFieldLoop = "source/Simulation/SimulationMain/magnetoHD/FieldLoop"
   SimMHDKH = "source/Simulation/SimulationMain/magnetoHD/KHmhd"
   SimMHDOT = "source/Simulation/SimulationMain/magnetoHD/OrszagTang"
   SimMHDMagRecon = "source/Simulation/SimulationMain/magnetoHD/MagneticRecon"
   SimMHDReconnection = "source/Simulation/SimulationMain/magnetoHD/Reconnection"
   SimMHDRotor = "source/Simulation/SimulationMain/magnetoHD/Rotor"
   SimMHDRT = "source/Simulation/SimulationMain/magnetoHD/RTmhd"
   SimNei = "source/Simulation/SimulationMain/NeiTest"
   SimNoh = "source/Simulation/SimulationMain/Noh"
   SimOrbit = "source/Simulation/SimulationMain/Orbit"
   SimPancake = "source/Simulation/SimulationMain/Pancake"
   SimPoisRefinement = "source/Simulation/SimulationMain/PoisTest_particleBasedRefine"
   SimProto = "source/Simulation/SimulationMain/ProtoPD"
   SimRT = "source/Simulation/SimulationMain/RTFlame"
   SimTFlame = "source/Simulation/SimulationMain/TurbFlame"
   SimRemap = "source/Simulation/SimulationMain/RemapGrid"
   SimSamTest = "source/Simulation/SimulationMain/SamraiTest"
   SimSBlast = "source/Simulation/SimulationMain/SBlast"
   SimSedov = "source/Simulation/SimulationMain/Sedov"
   SimSedovSG = "source/Simulation/SimulationMain/SedovSelfGravity"
   SimShock = "source/Simulation/SimulationMain/ShockCyl"
   SimSod = "source/Simulation/SimulationMain/Sod"
   SimSodRHD = "source/Simulation/SimulationMain/SodRHD"
   SimSodSphere = "source/Simulation/SimulationMain/SodSpherical"
   SimSodStep = "source/Simulation/SimulationMain/SodStep"
   SimSoundwave = "source/Simulation/SimulationMain/Soundwave"
   SimStir = "source/Simulation/SimulationMain/StirTurb"
   SimStirEIP = "source/Simulation/SimulationMain/StirTurbEIP"
   SimEosLite = "source/Simulation/SimulationMain/TestEosLite"
   Sim2Gamma = "source/Simulation/SimulationMain/TwoGamma"
   SimUTBurn = "source/Simulation/SimulationMain/unitTest/Burn"
   SimUTRayP = "source/Simulation/SimulationMain/unitTest/RayPath"
   SimUTPFFT = "source/Simulation/SimulationMain/unitTest/PFFT"
   SimUTPParab = "source/Simulation/SimulationMain/unitTest/PFFT_PoissonParabolae"
   SimUTPFFTXYpZnMC = "source/Simulation/SimulationMain/unitTest/PFFT_XYperZneu_MC"
   SimUTPFFTXYpnG2DMC = "source/Simulation/SimulationMain/unitTest/PFFT_XYperneu_GenDir2D_MC"
   SimUTUGReo = "source/Simulation/SimulationMain/unitTest/UGReordered"
   SimVortexRings = "source/Simulation/SimulationMain/VortexRings"
   SimWD = "source/Simulation/SimulationMain/WD_def"
   SimWDexp = "source/Simulation/SimulationMain/WD_def_exp"
   SimWDinert = "source/Simulation/SimulationMain/WD_def_inert"
   SimWDdet = "source/Simulation/SimulationMain/WD_det"
   SimWDjet = "source/Simulation/SimulationMain/WD_jet"
   SimWindTunnel = "source/Simulation/SimulationMain/WindTunnel"
   chomboUnitTests = "source/Grid/GridMain/Chombo/wrapper/unit_tests/"
 
   # dont check files in these directories for problems

   generalexcludes = [pmTwo, pmFour0Kernel]
   
   excludes = [pmTwo, pmFour0Kernel, pmFourDev, pmInterp,
       Samrai, SamraiIO, SPPM, SlowFlow,
       paraBurnOld, hgFlash2, hgFlash2fft,
       MultigridDiffHgADI, MultigridDiffHgFFT,
       EosNuc,
       cosmoKernel,
       ST_H,
       DiffUG, DiffCG,
       ImBound, SolidMechanics,
       IncompNSVarDens, MultigridFFT, MultigridExperimental, BIPCGStabVarcoeff,
       MultigridMC2,
       GridStructures,
       PartPredCorOld,
       SimAdvect, SimAdvectMass, SimBubLab, SimCellpb,
       SimCollShock, SimFlame3, SimFlameBub,
       SimDustCart,
       SimGalaxy, SimGalaxyM,
       SimClustSlos, SimClustSlosM, SimClustSlosMN,
       SimCoreColl,
       SimMHDAlfven, SimMHDBlastMHD,   SimMHDCloud,
       SimMHDDegen,  SimMHDKH,SimMHDMagRecon,SimMHDReconnection,
       SimMHDRT,
       SimINS, SimUTPFFTXYpZnMC, SimUTPFFTXYpnG2DMC,
       SimUTPParab,
       SimUTRayP,
       SimRemap,
       Sim2Gamma,
       SimFlash2Convert,
       SimGadget, SimLayer3, SimNoh, SimPoisRefinement, SimProto,
       SimTFlame,
       SimSamTest, SimShock, SimSoundwave,
       SimStirEIP, SimEosLite,
       SimVortexRings, SimWD, SimWDexp, SimWDinert,  SimWDdet, SimWDjet,
       SimNucOToRT, SimNuc2Grid,
       chomboUnitTests
       ]
                     
   def isInterfaceFile(self,lines):
       """
       Custom 'isComplex' function for some violation(s).
       Returns True if the file is too complex. Complex = is an interface module
       """
       if not lines:
          bfnint = self.regs["bfnInterface"]
          if bfnint.match(self.basename):
             return True
          else:
             return False
       begmod = self.regs["Mod"]
       endmod = self.regs["endMod"]
       begint = self.regs["endInterface"]
       inmod = 0
       for line in lines:
           if begmod.match(line):
              inmod = 1
           elif inmod and endmod.match(line):
              inmod= 0
           elif inmod and begint.match(line):
              return True
       return False

   # violations we currently handle
   violations = [ { "nickname" : "devComments", 
                    "fullname" : "DEV Comments",
                    "desc"     : "flag developer comments in source code"},
                  { "nickname" : "implicitNone", 
                    "fullname" : "IMPLICIT NONE",
                    "isComplex" : isInterfaceFile,
                    "desc"     : "Ensure that all subroutines have an IMPLICIT NONE statement",
                    "exceptions": ["^.*_.*[Ii]nterface[.]F90$"] },
                  { "nickname" : "intent",
                    "fullname" : "explicit declaration of parameter INTENT",
                    "desc"     : "The INTENT (IN,OUT,INOUT) must be explicity declared"},
                  { "nickname" : "save",
                    "fullname" : "Lonely SAVE",
                    "desc"     : "Using save by itself without specifying variables to save"},
                  { "nickname" : "noTabs",
                    "fullname" : "No TABS",
                    "desc"     : "Using TABS in Fortran source makes code non-portable. Consider replacing with ACHAR(9) "},
                  { "nickname" : "cleanSource",
                    "fullname" : "Source code which is 8bit clean and no control chars",
                    "desc"     : "Allow only ascii chars between 32 and 126 (and \\r\\n)"},
                  { "nickname" : "commonBlock",
                    "fullname" : "common Blocks are bad",
                    "desc"     : "Avoid using common Blocks. Exceptions must be hardcoded in script"},
                  { "nickname" : "fileName",
                    "fullname" : "sensible filenames",
                    "desc"     : 'name of source file must be name of the "main" subroutine or module in it'},
                  { "nickname" : "useOnly", 
                    "fullname" : "use USE always with ONLY (exceptions: use of some interface modules, use in _init files)",
                    "desc"     : "Use of USE (for data) without ONLY qualifier makes things harder to debug",
                    "exceptions": ["^.*_init[.]F90$"]},
                  { "nickname" : "FnArgsCheck", 
                    "fullname" : "Check all decl of functions are consistent",
                    "desc"     : "Ensure that all declarations of a function have same arguments (same names and order)"},
                  { "nickname" : "RoboDocCheck", 
                    "fullname" : "RoboDoc check",
                    "desc"     : "Checks for inconsistent robodoc header"},
                  { "nickname" : "deallocateCheck",
                    "fullname" : "balanced allocation / deallocation",
                    "desc"     : "Checks that anything for which memory is explicitly allocated is also deallocated."},
                  { "nickname" : "interfaces",
                    "fullname" : "consistency of API-level function calls / listings in API-level interface",
                    "desc"     : "Calls to API-level subroutines must be present in appropriate interface"}
                ]

   fixes = [ {"nickname" : "addImplNone",
              "fullname" : "Add implicit Nones",
              "isComplex": isInterfaceFile,
              "desc"     : "Add implicit nones for F90 subroutines which don't have them"},
             { "nickname": "TabsToSpace",
               "fullname": "Convert TABS to spaces",
               "desc"    : "Replace TABS in source with spaces"},
             { "nickname": "inclRoboDoc",
               "fullname": "Remove #includes in RoboDoc Headers",
               "desc"    : "Remove #includes in RoboDoc Headers"},
             { "nickname": "GenRoboDoc", 
               "fullname": "Generate RoboDoc header",
               "desc"    : "Generate Robodoc header if no header found"},
           ]

   # dictionary of information we already know based on name of argument
   # name of argument -> dictionary of info. When generating ROBODOC headers
   # this information in consulted to generate the description.
   # when info is not available the info for "none" is used
   RBinfo = { "none"     : "",
              "mype"     : "my Processor Number",
              "numprocs" : "total number of processors",
              "blockid"  : "ID of block in current processor",
              "blklimits": "array holding upper and lower index limits of interior block cells (no GC)",
              "blklimitsgc": "array holding the upper and lower index limits of an entire block (including GC)",
            }

   # list of directories we dont check or fix (we dont maintain code or too complex)
   nohandledirs = [ ]

   # types of file extensions we handle
   extensions = [ "C","F","F90" ]

   # commonly used regexps
   regexps = {
       "commonBlock" : r"^\s*\bcommon\s*[/](?P<name>.*)[/]",
       "subOrMod"    : r"^\s*\b((?:recursive\s+)?subroutine|module|program|\S*\s*function)\s+(?P<name>[-\w]*)",
       "Mod"         : r"^\s*\b(module)\s+(?P<name>[-\w]*)",
       "cFunction"   : r'^\s*\b(extern\s+"C"\s+)?\s*(?:\S+)?\s+(?:FTOC[(])?(?P<name>[-\w]*)[)]?',
       "devComments" : r"!\s*DEV.*$",
       "onlySave"    : r"^\s*\bsave\s*([!].*)?$",
       "useLine"     : r"^\s*\buse\s",
       # The next one matches a "use" statement followed by either a name ending in _interface, or having an ONLY.
       # Note that the next line allows cheating by putting a '!' right before the ",ONLY".
       "useOnlyLine" : r"^\s*\buse\s+\b([A-Za-z][A-Za-z0-9_]*\s*!?,\s*ONLY\s*:|[A-Z][A-Za-z0-9]*_interface\s*$)",
       "contains"    : r"^\s*contains\s*$",
       "begSub"      : r"^\s*\b(?:recursive\s+)?(?:subroutine|(?P<type>\S+)\s+function)\s+(?P<subname>[-\w]*)\s*\(",
       "fullSub"     : r"^\s*\b(?:recursive\s+)?(?:subroutine|\S+\s+function)\s+(?P<subname>[\w]*)\s*[(](?P<arglist>[^)]*)[)].*$", 
       # to completely parse subroutine decl call parseSubDecl
       "endSub"      : r"^\s*\bend\s+(?:subroutine|function)",
       "endMod"      : r"^\s*\bend\s+(?:module)",
       "subCall"     : r"(?:^|.*?\s+)call\s+(?P<subname>(?P<unitname>[a-zA-Z0-9]*?)(?:_\w+)?)\s*\(",
       "begInterface": r"^\s*\binterface(\s+(?P<name>\w+))?\s*$", 
       "begGenInterface": r"^\s*\binterface\s+(?P<name>\w+)\s*$", 
       "endInterface": r"^\s*\bend\s+interface",
       "bfnInterface": r"^[A-Za-z][A-Za-z0-9_]+_(interface|[a-z][a-z0-9]{1,12}[Ii]nterface)s?$",
       "intentLine"  : r"^[^!]*\b((INTENT\s*[(].*[)])|(pointer)).*$",
       "fullintLine" : r"^\s*(?P<typedecl>[^!]*\b((INTENT\s*[(].*[)])|(pointer)).*)\s*::(?P<ilist>[^!]*)(?:!.*)?$",
       "externalLine": r"^\s*external\s*(::)?\s*(?P<ilist>[^!]*)",
       "implNone"    : r"^\s*\bimplicit\s+none",
       "blankLine"   : r"^\s*(!.*)?$",
       "contLine"    : r"^\s*(&)?(?P<code>.*)[&]\s*(!.*)?$",
       "lastCLine"   : r"^\s*(&)?(?P<code>.*)\s*(!.*)?$",
       "pullArgs"    : r"\s*([^,(]+)(?:\([^()]*(?:\([^()]*\)[^()]*)*\))?\s*(?:,)?",
       "roboDoc"     : r"^!![*]+[a-z]+[*]+\s+\b(?P<name>.*)$",
       "rbSynArgs"   : r"::\s*([^\s,()]*)",
       "rbSynName"   : r".*?\b([^\s(]*)\s*[(]",
       "rbSynModName": r".*?\buse\s*([^\s(]*)\b",
       "rbInclude"   : r"^[!]{,2}\s*#(?:(?:include)|(?:inlcude))\s",
       "rbPullArgs"  : r"^\s*\b(?P<ilist>(?:(?:[^-\s:(])|(?:,\s*)|(?: *[(][^)]*[)]))*)\s*[-:]",
       "allocate"    : r".*?allocate\s*[(]\s*(.*?)\s*[(]",
       "deallocate"  : r".*?deallocate\s*[(]\s*(.*?)\s*[)]",
       "doubleQuote" : r'(?<!\\)".*?(?<!\\)"',
       "singleQuote" : r"(?<!\\)'.*?(?<!\\)'",
       "bfnData"     : r"^[A-Za-z][A-Za-z0-9_]+_(data|[a-z][a-z0-9]{1,2}[Dd]ata)$",
       "fnApi"       : r"^[A-Z]",
       "comment"     : r"^!.*$"

   }
   # pullArgs is an utterly compilcated regexp what is matches is the following
   # <whitespace><identifier><dimension>?<whitespace><comma>?
   # where the <comma> is optional, whitespace can be multiple whitespace chars (incl 0)
   # the dimension part is optional, and it will match any string which is delimited by ()
   # as long as the string between the () does not contain () nesting of more than one level
   #
   # i.e. <dimension> will match
   #  (5,4) and (ilogc:ihigc,jlogc:jhigc) and even 
   #  (extent(1,IAXIS):extent(2,JAXIS), extent(1,JAXIS):extent(2,JAXIS)) but not
   #  (foo(bar(1)):foo(bar(2)))
   #
   #  Writing a regexp to handle arbitrary levels of nesting is IMPOSSIBLE, you need
   #  Context Free Languages for that

   regs = {} # will be populated with compiled versions of regexps with re.I flag

   # fields to expect in a robodoc header
   RBfields = ["NAME","SYNOPSIS","DESCRIPTION","EXAMPLE","ARGUMENTS","NOTES",
               "PARAMETERS", "RESULT", "SIDE EFFECTS","SEE ALSO","AUTOGENROBODOC"]

   def __init__(self,flashHome,rptname):
       dict.__init__(self)
       if rptname == "-":
          self.rptname = "/dev/stdout"
       else: self.rptname = os.path.abspath(rptname)
       self.flashHome = flashHome
       excludes = [ os.path.abspath(os.path.join(flashHome,x)) for x in self.excludes ]
       self.excludes[:] = excludes
       self.methDict = {}
       self.filename = "" # the current filename we are working on
       self.complex_files = []
       self.complex_files_byMethod = []
       self.populate_check_methods()
       self.populate_fix_methods()
       self.compile_regexps()
       self.interfacesDict = InterfacesDict(flashHome)

   def compile_regexps(self):
       self.regs = {}
       for k,v in self.regexps.items():
           self.regs[k] = re.compile(v,re.I)

   def contractContLines(self,lines):
       """Given a list of lines returns a pair (# lines, line)
          line = is the first F90 statement line (removes continuation chars and comments)
          # lines is the number of lines from given list consumed to construct line
        """
       retlines = [] 
       rv = 0
       reg = self.regs["contLine"]
       m = reg.match(lines[0])
       while m:
           retlines.append(m.group("code"))
           rv += 1
           m = reg.match(lines[rv])
       else:
           m = self.regs["lastCLine"].match(lines[rv])
           retlines.append(m.group("code"))
       ans = " ".join(retlines)
       line = re.sub(r"[\n\r]","",ans)
       return (rv+1,line)
       
   def parseSub(self,lines):
       "Parses a subroutine declaration. Returns # of lines used for decl"
       data = {}
       line = ""
       m = self.regs["begSub"].match(lines[0])
       if not m: return ({},1) # not a sub
       # if you matched an "end function fubar" then ignore
       if m.group("type") == "end": return ({},1)
       # if no ( in lines[0] then the declaration is already over
       if lines[0].find("(") < 0:
          data["subname"] = m.group("subname")
          data["arglist"] = []
          return (data,1)
       (x,bigline) = self.contractContLines(lines)
       if bigline.find(")") < 0: # No ) found
          print "Some big parsing error while parsing %s" % self.filename
          sys.exit(1)
       fs = self.regs["fullSub"]
       m = fs.match(bigline)
       if not m: 
          return ({},1) # this is not subroutine at all
       data["subname"] = m.group("subname")
       data["arglist"] = filter(None, [y.strip() for y in 
                                       m.group("arglist").split(",")])
       return (data,x)
       
   def parseRB(self,lines):
       "Parses the ROBODOC header of a file"
       currfield = None
       ans = {}
       rv = {"syn_name": "","syn_args":[], "args": [], "problems": []}
       for key in self.RBfields: ans[key] = []
       for line in lines[1:]:
           if line[:2] != "!!": 
              rv["problems"].append("Found RoboDoc line not starting with !!")
              while line[:2] != "!!": line = "!"+line
           if line[:5] == "!!***": break
           line = line[2:].strip()
           if self.regs["rbInclude"].match(line) and currfield=="SYNOPSIS": # #include line
              rv["problems"].append("#include in Synopsis (RoboDoc header)")
              continue
           first = " ".join(re.split('\W+',line))
           if first in self.RBfields:
              currfield = first
              if currfield == "AUTOGENROBODOC":
                 rv["problems"].append("fix robodoc header descriptions and remove the AUTOGENROBODOC field")
              continue
           if not currfield: continue # no key to put it in
           line = line.replace("\n"," ")
           if line: 
              ans[currfield].append(line)
       synopsis = " ".join(ans["SYNOPSIS"])
       m = self.regs["rbSynName"].match(synopsis)
       if m:
          parts = {}
          for k in self.regs["rbSynArgs"].findall(synopsis): parts[k] = 1
          name = m.groups()[0]
          rv.update({"syn_name": name,"syn_args":parts.keys(), "args": []})
       else: # may be it is a module
          m = self.regs["rbSynModName"].match(synopsis)
          if m:
             name = m.groups()[0]
             rv.update({"syn_name": name, "syn_args": [],"args":[]})
       for line in ans["ARGUMENTS"]:
           m = self.regs["rbPullArgs"].match(line)
           if m: 
              # remove balanced parenthesis contents, i.e. dimension
              ilist = re.sub("[(][^)]*[)]","",m.group("ilist"))
              rv["args"].extend( [ x.strip() for x in ilist.split(",") ])
       args = {}
       for x in rv["args"]: args[x] = 1
       rv["args"][:] = args.keys()
       # do we have an arguments section but did not find any?
       if not rv["args"] and rv["syn_args"]:
          if " ".join(ans["ARGUMENTS"]).strip(): # found non-trivial arguments
             rv["problems"].append("Found argument section with no arguments")
       return rv

   def isComplex(self,lines):
       """
       Default 'isComplex' method for violations that don't name their own.
       Returns True if the file is too complex. Complex = has generic interface line
       """
       begint = self.regs["begGenInterface"]
       for line in lines:
           if begint.match(line):
              return True
       return False

   def print_complex(self,fd):
       fd.write("---------------------\n")
       fd.write("Skipped Complex Files\n")
       fd.write("---------------------\n")
       fd.write("\n")
       fd.write("   Files too complex and hence not checked/fixed.\n")
       fd.write("   Complex: Contains generic interface definitions\n")
       fd.write("\n\n")
       for f in self.complex_files: fd.write("   %s\n" % f)
       if (self.complex_files_byMethod):
          fd.write("\n")
          fd.write("   Files too complex and hence skipped for some checks/fixes.\n")
          fd.write("   Complex: check- or fix-specific criteria\n")
          fd.write("\n")
          for f in self.complex_files_byMethod: fd.write("   %s (skip %s %s)\n" % (f[1],f[0],f[2]) )
       fd.write("\n\n")

   def print_interfaceErrors(self, fd):
     if self.interfacesDict.errorsFound:
       fd.write("--------------------\n")
       fd.write("Errors in interfaces\n")
       fd.write("--------------------\n")
       fd.write("\n")

       missingInterfaces = [item for item in self.interfacesDict.keys() if
                            self.interfacesDict[item] == INTERFACE_FILE_MISSING]
       fd.write("The following API-Level directories are missing interface.F90 files:\n")
       fd.write("--------------------------------------------------------------------\n")
       fd.write("\n".join(missingInterfaces))
       fd.write("\n\n")

       for pathToAPILevelDir in self.interfacesDict.keys():
         subnamesDict             = self.interfacesDict[pathToAPILevelDir]
         missingDeclarations      = []
         caseMismatchDeclarations = []
         extraDeclarations        = []
         if isinstance(subnamesDict, dict):  # will be of type int if interface file was missing
           missingDeclarations      = [item for item in subnamesDict.keys() if
                                       subnamesDict[item] == INTERFACE_DECL_MISSING]
           caseMismatchDeclarations = [item for item in subnamesDict.keys() if
                                       subnamesDict[item] == INTERFACE_CASE_MISMATCH]
           extraDeclarations        = [item for item in subnamesDict.keys() if
                                       subnamesDict[item] == INTERFACE_DECL_EXTRA]

         pathToInterfaceFile = os.path.join(pathToAPILevelDir, (os.path.basename(pathToAPILevelDir) + "_interface.F90"))
         if missingDeclarations:
           fd.write("The following API-Level stub-files\n")
           fd.write("are not represented in \"%s\":\n" % pathToInterfaceFile)
           fd.write("--------------------------" + "-"*len(pathToInterfaceFile) + "\n")
           fd.write("\n".join(missingDeclarations))
           fd.write("\n\n")
         if caseMismatchDeclarations:
           fd.write("The capitalization patterns of the following interfaces from \"%s\"\n" % pathToInterfaceFile)
           fd.write("are not in agreement with that of the correspoding API-level stub-files:\n")
           fd.write("---------------------------------------------------------------------------------\n")
           fd.write("\n".join(caseMismatchDeclarations))
           fd.write("\n\n")
         if extraDeclarations:                    
           fd.write("The following subroutine interfaces from \"%s\"\n" % pathToInterfaceFile)
           fd.write("do not correspond to any API-level stub-file:\n")
           fd.write("-------------------------------------------" + "-"*len(pathToInterfaceFile) + "\n")
           fd.write("\n".join(extraDeclarations))
           fd.write("\n\n")

 
   def populate_fix_methods(self):
       self.fix_methods = {}
       for fix in self.fixes: 
           v = "fix_%(nickname)s" % fix
           # find default method to find those violations
           default = self.__class__.__dict__.get(v,None) # get the method whose name is x
           if not default: default = self.__class__.__dict__["default_fix"]
           for e in self.extensions:
               key = "%s_%s" % (v,e)
               self.fix_methods[key] = self.__class__.__dict__.get(key,default)

   def populate_check_methods(self):
       self.methods = {}
       for vio in self.violations: 
           v = vio["nickname"]
           # initialize entry to store violations of type v
           self[v] = [] # list of pairs (filename, violation_datastructure)
           # find method to print info regarding such violations
           pm = self.__class__.__dict__["default_print"] # default print method
           self.methods["%s_print_all" % v] = self.__class__.__dict__.get("%s_print_all" % v,None)
           self.methods["%s_print" % v] = self.__class__.__dict__.get("%s_print" % v,pm)
           # find default method to find those violations
           default = self.__class__.__dict__.get(v,None) # get the method whose name is x
           if not default: default = self.__class__.__dict__["default_violation"]
           for e in self.extensions:
               key = "%s_%s" % (v,e)
               self.methods[key] = self.__class__.__dict__.get(key,None)
               if not self.methods[key]: self.methods[key] = default
       
   def process(self,dirOrFileName,mode,methlist):
       pm = self.__class__.__dict__.get("%s_dir" % mode,None)
       if not pm: return None
       pm(self,os.path.expanduser(dirOrFileName),methlist)
       return 1

   def list_dir(self,dirOrFileName,methlist):
       print "\nNAMEOFMETHOD: Short description of method"
       for (k,dicts) in [("Check Methods",self.violations),
                         ("Fix Methods",self.fixes)]:
           print "\n%s\n%s" % (k,"-"*len(k))
           for d in dicts:
               print "\n%(nickname)s: %(fullname)s" % d
       print 

   def select_methods(self,dictmethods,methlist):
       if not methlist: return # all default state
       noValidMethods = True
       for d in dictmethods:
           if d["nickname"] not in methlist:
              d["disable"] = True
              noValidMethods = False
       if noValidMethods:
          print "No valid methods found. Methods given: %s" % ", ".join(methlist)
          print "Run with --mode=list to find valid methods"
          sys.exit(1)

   def check_file(self,filename):
       """Checks for all violations for the given file"""
       parts = string.split(os.path.basename(filename),".")
       self.basename = parts[0]
       if len(parts) > 1: self.ext = parts[1].upper()
       fd = open(filename)
       lines = fd.readlines()
       fd.close()
       for vio in self.violations:
           if vio.has_key("isComplex"):
              if vio["isComplex"](self,lines):
                 self.complex_files_byMethod.append( (vio["nickname"],filename,"check") )
                 continue
           elif self.isComplex(lines):
              self.complex_files.append(filename)
              return
           # else
           if vio.get("disable"): continue # dont do this violation
           v = vio["nickname"]
           key = v+"_"+self.ext
           check = True
           for regexp in vio.get("exceptions",[]): # for all prefixes in exceptions for this violation
               if re.match(regexp,filename): check = False
           if check:
              try:
                 ans = self.methods[key](self,lines)
              except:
                 print "Error while processing '%s' for file '%s'" % (key,filename)
                 raise
              if ans: self[v].append((filename,ans))

   def check_dir(self,dirOrFileName,methlist):
       self.select_methods(self.violations,methlist)
       self.process_dir(dirOrFileName,self.check_file)
       self.print_report()

   def edit_dir(self,dirOrFileName,methlist):
       self.select_methods(self.violations,methlist)
       self.process_dir(dirOrFileName,self.check_file)
       self.edit_violations(dirOrFileName)

   def fix_file(self,filename):
       """Checks for all violations for the given file"""
       parts = string.split(os.path.basename(filename),".")
       self.basename = parts[0]
       if len(parts) > 1: self.ext = parts[1].upper()
       fd = open(filename)
       lines = fd.readlines()
       fd.close()
       for vio in self.fixes:
           if vio.has_key("isComplex"):
              if vio["isComplex"](self,lines):
                 self.complex_files_byMethod.append( (vio["nickname"],filename,"fix") )
                 continue
           elif self.isComplex(lines):
              self.complex_files.append(filename)
              return
           # else
           if vio.get("disable"): continue # dont do this violation
           v = vio["nickname"]
           key = "fix_%s_%s" % (v,self.ext)
           check = True
           for pre in vio.get("exceptions",[]): # for all prefixes in exceptions for this violation
               if filename[:len(pre)] == pre: check = False
           if check:
              ans = self.fix_methods[key](self,lines)
              if not ans: continue # do the next fix
              # re write the file
              self.fix_report[v].append(filename)
              fd = open(filename,"w")
              fd.writelines(lines)
              fd.close()

   def fix_dir(self,dirOrFileName,methlist):
       self.select_methods(self.fixes,methlist)
       self.fix_report = {} # map fixname to list of files fixed
       for f in self.fixes:
           self.fix_report[f["nickname"]] = []
       self.process_dir(dirOrFileName,self.fix_file)
       fd = open(self.rptname,"w")
       self.print_complex(fd)
       for f in self.fixes:
           if f.get("disable"): continue
           p = len(f["fullname"])
           fd.write("-"*p)
           fd.write("\n%s\n" % f["fullname"])
           fd.write("-"*p)
           fd.write("\n\n%s\n\n" % f["desc"])
           fd.write("Files fixed:\n")
           ref = self.fix_report[f["nickname"]]
           if ref:
              for fn in ref:
                  fd.write("   %s\n" % fn)
           else: fd.write("   No files were fixed\n")
           fd.write("\n\n")
       fd.close()

   def process_dir(self,dirOrFileName,bound_method):

       def vfunc(flist,dname,fnames):
           # ignore directories starting with "."
           if dname[0] == ".": return
           fnames.sort()
           dontDescend = []
           for x in fnames:
               if x[0] == ".": # names starting with . (usually subdirs) to be ignored
                  dontDescend.append(x)
                  continue
               jname = os.path.join(dname,x)
               aname = os.path.abspath(jname)
               # ignore non -files
               if os.path.isdir(jname):
                  if aname in self.excludes: # ignore this sub-tree
                     dontDescend.append(x)
                     continue
               if os.path.isfile(jname):
                  ignore = 0
                  for prefix in self.excludes: # check if bad prefix
                      if aname[:len(prefix)] == prefix:
                         ignore = 1
                  if ignore == 1: continue
                  parts = string.split(x,".")
                  if len(parts) < 2: continue
                  if parts[1].upper() in self.extensions:
                     flist.append(jname)
           for x in dontDescend:
              fnames.remove(x)

       pwd = os.getcwd()
       bn = os.path.basename(dirOrFileName)
       os.chdir(os.path.dirname(dirOrFileName))

       if os.path.isdir(dirOrFileName):
          flist = []
          os.path.walk(bn,vfunc,flist)
          # now we have the list of all files to process
       else:
          # We assume 'dirOrFileName' refers to a file
          # Note that 'self.excludes' and 'self.extensions' are both
          # ignored in this case. That is, we process the explicitly-
          # named file regardless of other rules.
          flist = [dirOrFileName]
       for fname in flist: 
           self.filename = fname
           bound_method(fname)
       os.chdir(pwd)

   def print_report(self):
       fd = open(self.rptname,"w")
       self.print_complex(fd)
       self.print_interfaceErrors(fd)
       for vio in self.violations:
           if vio.get("disable"): continue
           v = vio["nickname"]
           fn = vio["fullname"]
           fd.write("-"*len(fn)+"\n")
           fd.write(fn+"\n")
           fd.write("-"*len(fn)+"\n")
           fd.write("\n")
           fd.write("Description: ")
           first = 1
           for line in vio.get("desc",fn).split("\n"):
               if not first: fd.write(" "*13)
               first = None
               fd.write(line+"\n")
           fd.write("\n\nViolations:\n")
           pm = self.methods["%s_print_all" % v]
           if pm: # we have a global print method
              pm(self,fd,"   ")
           else:
              pm = self.methods["%s_print" % v]
              for (fname,info) in self[v]:
                  fd.write("   %s\n" % fname)
                  pm(self,fd,"      ",info)
                  fd.write("\n")
           fd.write("\n")
       fd.close()

   def edit_violations(self,dirOrFileName):
      if os.path.isdir(dirOrFileName):
         dirname = dirOrFileName
      else:
         # we assume that 'dirOrFileName' refers to a file
         dirname = os.path.dirname(dirOrFileName)

      for vio in self.violations:
         if vio.get("disable"): continue
         v = vio["nickname"]
         fn = vio["fullname"]
         print "-"*len(fn)
         print fn
         print "-"*len(fn)
         print ""
         print "Description: "
         first = True
         for line in vio.get("desc",fn).split("\n"):
            if first:
               print line
               first = False
            else:
               print (" "*13) + line
         print ""
         print "Violations"
         for (fname,info) in self[v]:
            print fname
            for a in info:
               print "     " + str(a)
            if "__iter__" in dir(info[0]):
               try:
                  lineNo = int(info[0][0])
               except:
                  lineNo = None
            else:
               try:
                  lineNo = int(info[0])
               except:
                  lineNo = None
            ans = raw_input("Edit this file? [y/n]")
            ans = ans.replace("\n","")
            if ans.lower().startswith("y") or ans.lower == "":
               if lineNo:
                  cmd = ["emacs", "-nw", "+%s" % lineNo, os.path.join(os.path.dirname(dirname), fname)]
               else:
                  cmd = ["emacs", "-nw", os.path.join(os.path.dirname(dirname), fname)]
               try:
                  status = subprocess.check_call(cmd)
                  if status != 0:
                     print "There were some errors from", cmd
               except NameError:
                  p = popen2.Popen3(cmd)
                  p.wait()
               print


   ###################### methods for identifying violations and printing the report
   def default_violation(self,lines):
       return []

   def default_print(self,fd,prefix,info):
       fd.write(prefix+str(info)+"\n")
       
   def RoboDocCheck_F90(self,lines):
       rv = []
       if not lines:
          return ["Empty file, no ROBODOC header"]
       m = self.regs["roboDoc"].match(lines[0])
       if not m:
          return ["No ROBODOC header"]
       ans = self.parseRB(lines)
       if ans.has_key("problems"): rv.extend(ans["problems"])
       if not ans["syn_name"].startswith(self.basename):
          rv.append("Synopsis name is '%s', filename is '%s'" % (ans["syn_name"],self.basename))
          if not ans["syn_name"]: 
             rv.append("Did you forget to use `use' before module name in synopsis?")
       syn = ans["syn_args"]
       args = ans["args"]
       common = [x for x in syn if x in args]
       for x in common:
           syn.remove(x)
           args.remove(x)
       if syn:
          rv.append("Arguments found only in SYNOPSIS: %s" % ", ".join(syn))
       if args:
          rv.append("Arguments found only in ARGUMENTS: %s" % ", ".join(args))
       return rv
       
   def RoboDocCheck_print(self,fd,prefix,info):
       for a in info:
           fd.write("%s %s\n" % (prefix,a))

   def stripQuotesAndComments(self, line):
     line = line.strip()
     line = self.regs["doubleQuote"].sub("", line)
     line = self.regs["singleQuote"].sub("", line)
     line = self.regs["comment"].sub("", line)
     line = line.strip()
     return line

   def deallocateCheck_F90(self, lines):
     allocates  = []
     returnList = []
     counter = 0
     for line in lines:
       counter += 1
       line = self.stripQuotesAndComments(line)
       if line:
         m = self.regs["allocate"].search(line)
         if m:
           allocates.append((m.group(1), counter))
           continue
         # else
         m = self.regs["deallocate"].search(line)
         if m:
           sought = m.group(1)
           for allocate, lineNo in allocates[:]:
             if allocate == sought:
               allocates.remove((allocate, lineNo))

     if allocates:
       returnList = [lineNo for allocate, lineNo in allocates]

     return returnList


   def interfaces_F90(self, lines):
     subCall = self.regs["subCall"]
     ans = []
     i = 0
     while i < len(lines):
       line = self.stripQuotesAndComments(lines[i])
       if line:
         m = subCall.match(line)
         if m:
           unitname = m.group("unitname")
           subname  = m.group("subname")
           pathToAPILevelDir, subnamesDict = self.interfacesDict.getPathAndSubroutinesForUnit(unitname)
           if pathToAPILevelDir:
             correctlyCasedUnitname = self.interfacesDict.getCorrectlyCasedUnitname(unitname)
             pathToInterfaceFile    = os.path.join(pathToAPILevelDir, (correctlyCasedUnitname + "_interface.F90"))

             if subname not in subnamesDict.keys():
               correctlyCasedSubname = self.interfacesDict.getCorrectlyCasedSubname(subname)
               if correctlyCasedSubname:
                 # Something has been typed with poor attention to upper/lower case
                 ans.append((str(i), subname, (pathToInterfaceFile, correctlyCasedSubname), INTERFACE_CASE_MISMATCH))
               else:
                 # This subroutine is not in the interface, but it might actually
                 # be an internal (non-API-level) subroutine that only matched
                 # 'subCall' because of that re-object's case-insensitivity (which
                 # was necessary for catching case mismatches in the above clause).
                 # Now, if 'unitname' is equal to'correctlyCasedUnitname', we'll
                 # know this is definitely supposed to be an API-level subroutine
                 # call and that we've really found an error.
                 if unitname == correctlyCasedUnitname:
                   ans.append((str(i), subname, pathToInterfaceFile, INTERFACE_DECL_MISSING))
       i+=1
     return ans

   def interfaces_print(self, fd, prefix, info):
     for lineNo, subname, pathToInterfaceFile, error in info:
       if error == INTERFACE_CASE_MISMATCH:
         # unpack the tuple represented by 'pathToInterfaceFile' for this type of error
         pathToInterfaceFile, realInterfaceName = pathToInterfaceFile
         fd.write("line %s: The capitalization pattern of called subroutine \"%s\"\n" % (lineNo, subname))
         fd.write("does not match the corresponding subroutine in the interface file \"%s\",\n" % pathToInterfaceFile)
         fd.write("where it appears as \"%s\".\n" % realInterfaceName)
         fd.write("\n")
       elif error == INTERFACE_DECL_MISSING:
         fd.write("line %s: The called subroutine \"%s\" is missing from interface file \"%s\".\n" % (lineNo, subname, pathToInterfaceFile))


   def fileName(self,lines):

       def contract(name):
           return filter(lambda x: x != "_", name).lower() # remove _ from name

       if self.ext in ["C","c"]:
          reg = self.regs["cFunction"]
       else: reg = self.regs["subOrMod"]
       snamereg = re.compile("^%s.*$" % contract(self.basename))
       ok = None
       for line in lines:
           m = reg.match(line)
           if m: 
              cname = contract(m.group("name"))
              if snamereg.match(cname): ok = 1
       if not ok:
          return self.basename
       else: 
          return None

   def fileName_print(self,fd,prefix,info):
       fd.write("%s subroutine (or module or C function) starting with %s not found\n" % (prefix,info))

   def commonBlock_F90(self,lines):
       exceptions = ["source/Grid/GridMain/paramesh/interpolation/Paramesh3/prolong/amr_prolong_gen_work1_fun.F90",
                     "source/Grid/GridMain/paramesh/interpolation/Paramesh3/prolong/amr_prolong_gen_unk1_fun.F90",
                     "source/Grid/GridMain/paramesh/interpolation/Paramesh4/prolong/amr_prolong_gen_work1_fun.F90",
                     "source/Grid/GridMain/paramesh/interpolation/Paramesh4/prolong/amr_prolong_gen_unk1_fun.F90"]

       if self.filename in exceptions:
          return []

       reg = self.regs["commonBlock"]
       lno = 1
       ans = []
       for line in lines:
           m = reg.match(line) # is there a common block
           if m: ans.append(m.group("name"))
           lno += 1
       return ans

   def commonBlock_print(self,fd,prefix,info):
       fd.write("%s Common Blocks: %s\n" % (prefix,string.join(info,", ")))
       
   def devComments_F90(self,lines):
       reg = self.regs["devComments"]
       lno = 1
       ans = []
       for line in lines:
           m = reg.search(line) # found a DEV comment
           if m:
              ans.append((lno,line[m.start():m.end()]))
           lno += 1
       return ans

   def devComments_F(self,lines):
       return self.devComments_F90(lines)
  
   def devComments_print(self,fd,prefix,info):
       for (n,com) in info:
           fd.write("%sLine %3d : %s\n" % (prefix,n,com))

   def save_F90(self,lines):
       ans = []
       lno = 1
       reg = self.regs["onlySave"]
       for line in lines:
           if reg.match(line): # found a standalone save
              ans.append(lno)
           lno += 1
       return ans

   def save_print(self,fd,prefix,info):
       for n in info:
           fd.write("%sLine %3d\n" % (prefix,n))

   def cleanSource_F90(self,lines):
       ans = [] # return list of line,pos,bad char
       lno = 0
       for line in lines:
           lno += 1
           p = 1
           for c in line:
               if c in ["\n","\r"]: continue
               if (ord(c) < 32) or (ord(c) > 126):
                  ans.append((lno,p,c))
               p += 1
       return ans

   def cleanSource_print(self,fd,prefix,info):
       for l,p,c in info:
           fd.write("%sLine %3d Position %d Char %s(%d)\n" % (prefix,l,p,c,ord(c)))

   def noTabs_F90(self,lines):
       ans = []
       lno = 1
       for line in lines:
           p = line.find("\t")
           if p >= 0: 
              ans.append((lno,p))
           lno += 1
       return ans

   def noTabs_print(self,fd,prefix,info):
       for n,p in info:
           fd.write("%sLine %3d Position %d\n" % (prefix,n,p))

   def useOnly_F90(self,lines):
       """Search for lines starting with USE and not followed by ONLY:"""
       useline = self.regs["useLine"]
       useonlyline = self.regs["useOnlyLine"]
       badlnos = []
       lno = 0
       for line in lines:
           lno += 1
           if useline.match(line) and not useonlyline.match(line): # found a bad line
              badlnos.append(lno)
       return badlnos
              
   def useOnly_print(self,fd,prefix,info):
       for n in info:
           fd.write("%sLine %3d\n" % (prefix,n))

   def implicitNone_F90(self,lines):
       begsub = self.regs["begSub"]
       endsub = self.regs["endSub"]
       contains = self.regs["contains"]
       impnone = self.regs["implNone"]
       badsubnames = []
       currsub = None # current subroutine
       currcontainingsub = None
       emptysub = 1 # is the subroutine empty?
       foundimpnone = 0
       for line in lines:
           if not currsub:
              m = begsub.match(line)
              if m:
                 # found beginning of sub
                 if not currcontainingsub: foundimpnone = 0
                 emptysub = 1
                 currsub = m.group("subname")
              elif currcontainingsub:
                 m = endsub.match(line)
                 if m:
                    foundimpnone = 0
                    emptysub = 1
                    currcontainingsub = None
           else:
              m = impnone.match(line)
              if m:
                 foundimpnone = 1
                 emptysub = None
                 continue
              else:
                 m = endsub.match(line)
                 if not m: 
                    m = contains.match(line)
                    if m:               # found CONTAINS
                       if not currcontainingsub:
                          if not foundimpnone: # no impNone in a subroutine with a CONTAINS statement
                             badsubnames.append(currsub + " (CONTAINS others)")
                          currcontainingsub = currsub
                          currsub = None
                       foundimpnone = 1 # until outer level subroutine ends, maybe
                       continue
                    else:
                       if len(line)>4: # heuristic for non-trivial code (we dont check comments)
                          emptysub = None
                       continue
                 if not foundimpnone and not emptysub: # no impNone in a non-trivial subroutine
                    if not currcontainingsub:
                       badsubnames.append(currsub)
                 foundimpnone = 0
                 emptysub = 1
                 currsub = None
       return badsubnames

   def implicitNone_print(self,fd,prefix,info):
       for name in info:
           fd.write("%s Subroutine %s\n" % (prefix,name))

   def FnArgsCheck_F90(self,lines):
       ans = []
       resub  = self.regs["begSub"]
       if not self.methDict.has_key("FnArgsCheck"): 
          self.methDict["FnArgsCheck"] = {}
       mdict = self.methDict["FnArgsCheck"] # store info across calls in this dict
       length = len(lines)
       dic = {}
       ctr = 0
       while ctr < length:
           line = lines[ctr]
           ctr += 1
           m = resub.match(line)
           if m:
              (dic,ctrinc) = self.parseSub(lines[ctr-1:])
              ctr += ctrinc - 1
              if not dic: continue # if bogus match proceed
              name = dic["subname"]
              args = [x.strip().lower() for x in dic["arglist"]]
              # mdict[name] is a list of pairs (filename,decl)
              if not mdict.has_key(name): 
                 mdict[name] = []
              # the different declarations we have for subroutine name
              decs = [ x[1] for x in mdict[name] ] 
              if args in decs: continue
              # we have a new declaration
              mdict[name].append( (self.filename, args))
              if len(mdict[name]) > 1: ans = name
       return ans
       
   def FnArgsCheck_print_all(self,fd,prefix):
       if not self.methDict.has_key("FnArgsCheck"): 
          return
       mdict = self.methDict["FnArgsCheck"]
       for (k,v) in mdict.items():
           if len(v) > 1: # more than one decl  for a function
              fd.write("\n%sFunction %s has multiple declarations\n" % (prefix,k))
              for (fname,args) in v:
                  fd.write("%s   %s\n%s      %s\n" % (prefix,fname,prefix,", ".join(args)))

   def intent_F90(self,lines):
       begsub   = self.regs["begSub"]
       intline  = self.regs["intentLine"]
       fullint  = self.regs["fullintLine"]
       external = self.regs["externalLine"]
       endsub   = self.regs["endSub"]
       pullargs = self.regs["pullArgs"]
       currsub = None # current subroutine
       noints = []
       dic = {}
       length = len(lines)
       ctr = 0
       ans = [] # list of pairs of subname, variables without INTENT
       while ctr < length:
           line = lines[ctr]
           ctr += 1
           if not currsub: # searching for start of subroutine
              m = begsub.match(line)
              if not m: continue
              (dic,ctrinc) = self.parseSub(lines[ctr-1:])
              ctr += (ctrinc - 1)
              if not dic: continue
              currsub = dic["subname"]
              noints = dic["arglist"]
           else: # we are inside a subroutine
              m = intline.match(line) # an INTENT line?
              if m: # remove vars for which INTENT has been found
                 (ctrinc,bigline) = self.contractContLines(lines[ctr-1:])
                 ctr += ctrinc - 1
                 m = fullint.match(bigline)
                 if not m: # we have a fake intent/pointer line (may be inside a string)
                    ctr -= ctrinc - 1
                    break
                 ilist = m.group("ilist")
                 all = pullargs.findall(ilist)
                 for x in all:
                     try:
                        noints.remove(x.strip()) 
                     except: pass
              else: # variables marked "external" don't need to declare intent
                 m = external.match(line)
                 if m: # remove vars which are marked "external"
                    (ctrinc,bigline) = self.contractContLines(lines[ctr-1:])
                    ctr += ctrinc - 1
                    m = external.match(bigline)
                    if not m: # we have a fake external line (may be inside a string)
                       ctr -= ctrinc - 1
                       break
                    ilist = m.group("ilist")
                    all = pullargs.findall(ilist)
                    for x in all:
                        try:
                           noints.remove(x.strip())
                        except: pass
                 else: # are we in end of sub
                    m = endsub.match(line)
                    if m: # end of subroutine found
                       if noints: ans.append((currsub,noints))
                       currsub = None
                       noints = []
       return ans

   def intent_print(self,fd,prefix,info):
       for name,vars in info:
           v = ", ".join(vars)
           fd.write("%s Subroutine: %s \n%s Variables: %s\n" % (prefix,name,prefix,v))

   ################## methods for fixes 

   def default_fix(self,lines):
       return None

   def fix_TabsToSpace(self,lines):
       rv = None
       TABSIZE = 6 # place tabstops every TABSIZE chars
       for x in range(len(lines)):
           conv = lines[x].expandtabs(TABSIZE)
           if conv != lines[x]:
              lines[x] = conv
              rv = 1
       return rv

   ## bunch of source code from after declaration to end subroutine
   # insert implicit none in appropriate place
   def insertImplNone(self,lines,start):
       length = len(lines)
       ctr = start
       bline = self.regs["blankLine"]
       uline = self.regs["useLine"]
       cline = self.regs["contLine"]
       while ctr < length:
             line = lines[ctr]
             ctr += 1
             if bline.match(line): continue
             if uline.match(line): # found a use
                while cline.match(line): # it continues to next line
                      line = lines[ctr]
                      ctr += 1
                continue
             lines[ctr-1:ctr-1] = ['implicit none !! Added by fix script\n']
             break

   def fix_inclRoboDoc(self,lines):
       outlines = []
       rv = None
       currfield = None
       for ctr in range(len(lines)):
           line = lines[ctr]
           if line[:2] != "!!":  # end of RoboDoc reached
              outlines.extend(lines[ctr:])
              break
           outlines.append(line) # add it anyway for now
           line = line [2:]
           first = " ".join(re.split('\W+',line)).strip()
           if first in self.RBfields:
              currfield = first
              continue
           if self.regs["rbInclude"].match(line) and currfield=="SYNOPSIS": # found a fix to make
              del outlines[-1] # remove the line we just added
              rv = 1
       lines[:] = outlines
       return rv
              
   def isF90DataFile(self,lines):
       """
       Custom 'isComplex' function for some violation(s).
       Returns True if the file is too complex. Complex = is an interface module
       """
       bfndata = self.regs["bfnData"]
       namematch = bfndata.match(self.basename)
       if namematch:
          if not lines: return True
       else:
          if not lines: return False
       begmod = self.regs["Mod"]
       endmod = self.regs["endMod"]
       begsub = self.regs["fullSub"]
       callsub = self.regs["subCall"]
       begallo = self.regs["allocate"]
       begint = self.regs["begInterface"]
       inmod = 0
       hasmod = 0
       for line in lines:
          if begsub.match(line):
             return False
          if callsub.match(line):
             return False
          if begallo.match(line):
             return False
          if begmod.match(line):
             inmod = 1
             hasmod = 1
          elif inmod and endmod.match(line):
             inmod= 0
          elif begint.match(line):
             return False
       if namematch:
          return True
       elif hasmod:
          return True
       return False

   def fix_GenRoboDoc_F90(self,lines):
       # no fix needed if header found
       if lines:
          m = self.regs["roboDoc"].match(lines[0])
       else:
          m = ''
       if m: return None

       # use FnArgsCheck to check for subroutine declaration
       self.methDict["FnArgsCheck"] = {}
       if not lines:
          return self.fix_GenRoboDocAndSkelForEmpty_F90(lines)
       self.FnArgsCheck_F90(lines)
       # mdict is a dictionary mapping FnName to (filename,listofargs)
       mdict = self.methDict["FnArgsCheck"]
       # search for an entry for name equal to filename
       if not mdict.has_key(self.basename): # no subroutine with name = filename
          return None # dont fix this yet, this will be flagged by other checks

       # reparse code to find the types of the arguments
       # and populate rbinfo
       args = mdict[self.basename][0][1]
       rbinfo = {}
       for k,v in self.RBinfo.items(): rbinfo[k] = {"type": "", "desc":v }
       inmysub = None
       ctr = 0
       length = len(lines)
       intline = self.regs["intentLine"]
       fullint= self.regs["fullintLine"]
       begsub = self.regs["begSub"]
       endsub = self.regs["endSub"]
       pullargs = self.regs["pullArgs"]
       while ctr < length:
           line = lines[ctr]
           ctr += 1
           m = begsub.match(line)
           if m and m.group("subname") == self.basename: # found my sub
              inmysub = 1
           m = endsub.match(line)
           if m and inmysub: # found end of sub, quit loop
              inmysub = None 
              break
           if not inmysub: continue
           # scan for intent line and pull out declaration
           m = intline.match(line) # an INTENT line?
           if m: # remove vars for which INTENT has been found
              (ctrinc,bigline) = self.contractContLines(lines[ctr-1:])
              ctr += ctrinc - 1
              m = fullint.match(bigline)
              if not m: # we got a fake intent line
                 ctr -= ctrinc - 1
                 break
              ilist = m.group("ilist")
              typedecl = m.group("typedecl")
              for x in pullargs.findall(ilist):
                  x = x.strip().lower()
                  if not rbinfo.has_key(x): 
                     rbinfo[x] = {}
                     rbinfo[x].update(rbinfo["none"]) # copy the defaults
                  typedecl = re.sub(" {2,}"," ",typedecl)
                  typedecl = re.sub(", *intent *[(](?P<intentkw>(in|out|inout))[)] *","(\g<intentkw>)",typedecl)
                  rbinfo[x]["type"] = re.sub(" {2,}"," ",typedecl)

       # Generate the robodoc header
       rdoc = [] # before printing we will add !! in front of all lines
       rdoc.append("****if* %s" % self.filename[:self.filename.find(".F90")])
       rdoc.extend([""," NAME","","  %s" % self.basename,""," SYNOPSIS",""])
       argprefix = " "*(2+5+len(self.basename))
       synlist = [] 
       arglist = []
       for arg in args:
           ai = rbinfo.get(arg.lower(),rbinfo["none"])
           synlist.append(" %s%s :: %s," % (argprefix,ai["type"],arg))
           arglist.append("")
           arglist.append("   %s : %s" % (arg,ai["desc"]))
       if not synlist: synlist = [""] # always non-empty
       synlist[0] = "  call %s(%s" % (self.basename,synlist[0].lstrip())
       if synlist[-1][-1] == ",":
          synlist[-1] = "%s)" % synlist[-1][:-1]
       else: synlist[-1] += ")"
       if not arglist: arglist = ["","   No arguments"] 
       rdoc.extend(synlist)
       rdoc.extend([""," DESCRIPTION","",""," ARGUMENTS"])
       rdoc.extend(arglist)
       rdoc.extend([""," AUTOGENROBODOC",""])
       rdoc.extend(["","***"])
       lines.insert(0,"\n")
       rdoc.reverse()
       for line in rdoc:
           lines.insert(0,"!!"+line+"\n")
       return 1

   def fix_GenRoboDocAndSkelForEmpty_F90(self,lines):
       # Generate the robodoc header for a completely empty file
       rdoc = [] # before printing we will add !! in front of all lines
       afilename = os.path.abspath(self.filename)
       srcpos = afilename.find("source/")
       if srcpos >= 0:
          filename = afilename[srcpos:]
       else:
          filename = self.filename
       synlist = [""] # always non-empty
       if self.isInterfaceFile(lines):
          stmt = 'use'
          isApi = self.regs["fnApi"]
          if isApi.match(self.basename):
             desig = 'h'
          else:
             desig = 'ih'
          synlist[0] = "  %s %s" % (stmt,self.basename)
       elif self.isF90DataFile(lines):
          stmt = 'use'
          desig = 'if'
          synlist[0] = "  %s %s, ONLY: ..." % (stmt,self.basename)
       else:
          stmt = 'call'
          desig = 'if'
          synlist[0] = "  %s %s(%s" % (stmt,self.basename,synlist[0].lstrip())
          if synlist[-1][-1] == ",":
             synlist[-1] = "%s)" % synlist[-1][:-1]
          else: synlist[-1] += ")"
       rdoc.append("****%s* %s" % (desig, filename[:filename.find(".F90")]))
       rdoc.extend([""," NAME","","  %s" % self.basename,""," SYNOPSIS",""])
       rdoc.extend(synlist)
       rdoc.extend([""," DESCRIPTION","",""," ARGUMENTS"])
       arglist = ["","   "]
       rdoc.extend(arglist)
       rdoc.extend([""," AUTOGENROBODOC",""])
       rdoc.extend(["","***"])
       lines.insert(0,"\n")
       rdoc.reverse()
       for line in rdoc:
           lines.insert(0,"!!"+line+"\n")
       if stmt == 'call':
          lines.append("subroutine %s()\n" % self.basename)
          lines.append("  implicit none\n\n")
          lines.append("end subroutine %s\n" % self.basename)
       else:
          lines.append("Module %s\n" % self.basename)
          lines.append("  implicit none\n\n")
          if desig == 'if':
             lines.append("  integer, save ::\n")
             lines.append("  logical, save ::\n")
             lines.append("  real, save ::\n")
             lines.append("\n")
          else:
             lines.append("  interface\n")
             unt = self.basename
             uspos = unt.find("_")
             if uspos > 0:
                unt = unt[0:uspos+1]
             lines.append("     subroutine %s\n" % unt)
             lines.append("     end subroutine %s\n" % unt)
             lines.append("  end interface\n\n")
          lines.append("end Module %s\n" % self.basename)
       return 1

   def fix_addImplNone_F90(self,lines):
       begsub = self.regs["begSub"]
       endsub = self.regs["endSub"]
       impNone= self.regs["implNone"]
       foundimpNone = 0
       dic = {}
       outlines = []
       rv = None
       currsub = [] # lines for current subroutine
       length = len(lines)
       ctr = 0
       subdeclen = 0
       while ctr < length:
           line = lines[ctr]
           ctr += 1
           if begsub.match(line): # we have a new subroutine
              currsub = []
              (dic,subdeclen) = self.parseSub(lines[ctr-1:])
              currsub.extend(lines[ctr-1:ctr+subdeclen-1])
              ctr += subdeclen - 1
              foundimpNone = 0
           elif impNone.match(line):
              currsub.append(line)
              foundimpNone = 1
           elif endsub.match(line):
              currsub.append(line)
              if foundimpNone == 0:
                 self.insertImplNone(currsub,subdeclen)
                 rv = 1
              outlines.extend(currsub)
              currsub = []
           else:
              if currsub: # in middle of a routine
                 currsub.append(line)
              else: # inbetween routines
                 outlines.append(line)
       lines[:] = outlines
       return rv


########################### START SCRIPT ################################
def usage():
   print >> sys.stderr, "Usage: codecheck.py --mode=<mode> --report-file=<report-file> --target=<target> --method=<method>"
   print >> sys.stderr 
   print >> sys.stderr, " <mode> is one of 'check', 'fix', 'list', 'edit'. Default 'check'."
   print >> sys.stderr, "     List mode lists the methods which can be used."
   print >> sys.stderr 
   print >> sys.stderr, " <report-file> is where the report will be written."
   print >> sys.stderr, "     If not specified, it defaults to 'docs/designDocs/violations.txt' for check mode."
   print >> sys.stderr, "     - implies stdout"
   print >> sys.stderr 
   print >> sys.stderr, " <target> If a directory, all files underneath it will be checked (recursively)."
   print >> sys.stderr, "     If a file, only that single file will be checked."
   print >> sys.stderr, "     if not specified, defaults to FLASH_HOME/source"
   print >> sys.stderr 
   print >> sys.stderr, " <method> determines which checks/fixes should be applied."
   print >> sys.stderr, "     If not specified defaults to all methods."
   sys.exit(1)


def main():
   # codeCheck.py is FLASH_HOME/tools/scripts/codeCheck.py
   # FLASH_HOME is the parent of directory having setup.py 
   flashHome = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))
   dirOrFileName = os.path.join(flashHome,"source")
   # command line args
   opts = ["mode=","report-file=","target=","method="]
   # default report filename
   fname = os.path.join(flashHome,"docs","designDocs","violations.txt")
   # default mode: choices check,fix
   mode = "check"
   meths = []
   try:
      (optlist, args) = getopt.getopt(sys.argv[1:],None,opts)
   except:
      usage()
   # no arguments needed
   if args: usage() 
   for (k,v) in optlist:
       if k == "--mode":
          mode = v
       elif k == "--report-file":
          fname = v
       elif k == "--target":
          dirOrFileName = os.path.abspath(v)
       elif k == "--method":
          meths.extend(v.split(","))
   cc = checkCode(flashHome,fname)
   if not cc.process(dirOrFileName,mode,meths):
      print >> sys.stderr, "Invalid mode %s\n" % mode
      usage()

########################### END SCRIPT #######################################
if __name__=='__main__':
    main()

