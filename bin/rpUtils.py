
# code related to Runtime Parameter Stuff

__all__ = [ "RPInfo", 
          ]

##############################################

import globals
from globals import *   # GVars, opts and SetupError
from lazyFile import *  # for LazyFile class


import string, re, sys, UserDict, os.path

###############################################################

class RPInfo(UserDict.UserDict):
   # Class which keeps track of run time parameters and related info
   # To each runtime parameter we keep track of the following info
   #
   # "NAME" of the parameter
   # "TYPE" - REAL, INTEGER, BOOLEAN, STRING
   # "VALUE" - Initial value of parameter
   # "CONST" - Boolean saying if this value is constant
   # "COMMENT" - Info regarding parameter
   # "RANGE" - [ list of range specifications ]
   #           range spec = {"min":min_val,"max":max_val} or "STRING" 
   #
   # self is a dictionary mapping "LOCATION","NAME" to other info
   # self.locations is a dictionary mapping "NAME" to list of locations defining it
   # If warn is non-empty
   #    duplicate names are considered WARNING
   # else
   #    duplicate (NAME,LOCATION) is WARNING else it is OK
   #
   # If it is a warning then RP not added to database

   def __init__(self,warn=1):
       UserDict.UserDict.__init__(self)
       self.typemap = {'REAL': 'real',
                       'INTEGER': 'integer',
                       'STRING':  'character(len=MAX_STRING_LENGTH)',
                       'BOOLEAN':'logical',
                       'DOC': ""}
       self.warn = warn
       self.locations = {}
       self.rangeToregexp = re.compile(r"to|TO|[.][.][.]")
       self.numregexp = re.compile(r"\s*(?P<min>\S*)\s*(?:to|TO|[.][.][.])\s*(?P<max>\S*)\s*")
       self.numre= re.compile(r"^[-+]?[0-9]*([.]([0-9]*([eE][-+]?[0-9]+)?)?)?$")
       self.numExpoRe= re.compile(r"^[-+]?0*((?P<before>[0-9]+)[.]?|(?P<before2>[0-9]*)[.](?P<aft0>0*)[0-9]+)[eEdD](?P<expo>[-+]?[0-9]+)$")
       self.docIgnoreRe = re.compile(r"^\s*__IGNORE__\s*$")

   def parseRange(self,range,typ,name):

       def checkNum(what):
           if not what: return 1
           if what=="TINY": return 1
           if what=="-TINY": return 1
           # strip trailing +,-
           if what[-1] in ["+","-"]: what = what[-1]
           if self.numre.match(what):
              return 1
           else: return None

       rv = []
       if typ not in ["REAL","INTEGER","STRING"]:
          return rv
       if not range:
          return rv
       for rangespec in range.split(","):
          
          rangespec = rangespec.strip()
          if not rangespec: continue
          if (self.rangeToregexp.search(rangespec) >= 0) and (typ in ["REAL","INTEGER"]): # contains "to" as a substring
             spec = {"min":None,"max":None}
             m = self.numregexp.match(rangespec)
             spec.update(m.groupdict())
             rv.append(spec)
          else:
             if rangespec[0] == rangespec[-1] and rangespec[0] in ["'",'"']: 
                rangespec = rangespec[1:-1]
             if typ in ["REAL","INTEGER"]: # checking for membership in degenerate interval
                rv.append( {"min":rangespec,"max":rangespec})
             elif typ in ["STRING"]: 
                rv.append(rangespec)
       # check if given numbers are valid
       if typ in ["REAL","INTEGER"]:
          for rs in rv: # for each range specification
            if not checkNum(rs["min"]): # min is bad
               raise SetupError("minimum value '%s' is illformed in parameter '%s'" % (rs["min"],name))
            if not checkNum(rs["max"]): # min is bad
               raise SetupError("maximum value '%s' is illformed in parameter '%s'" % (rs["max"],name))
       return rv

   def addRP(self,name,type="",value="",const=None,location="",comment="",range=None):
       if self.warn:
          if self.locations.has_key(name):
             loc2 = self.locations[name][0]
             val2 = self[(loc2,name)]["VALUE"]
             typ2 = self[(loc2,name)]["TYPE"]
             if (val2 == value) and (typ2 == type):
                # already present ignore this call
                if (self[(loc2,name)]["RANGE"] == self.parseRange(range,type,name)): return
             GVars.out.push()
             quit = 1 # quit this call
             sP = globals.simulationPrefix
             if sP.endswith(os.sep): sP = sP[:-1]
             # check if either of them is in SimlulationUnit. If so that overrides (exclusive OR)
             if (loc2.startswith(sP) and not location.startswith(sP)) or \
                    (not loc2.startswith(sP) and location.startswith(sP)):
                GVars.out.put('\nINFO: Parameter %s defined in both\n%s (default %s) and \n%s (default %s)'
                        %(name, loc2, val2, location, value),globals.IMPINFO)
                GVars.out.put("Simulation instance overrides; removing other instance.",globals.IMPINFO)
                if location.startswith(sP): # remove existing details
                   del self[(loc2,name)]
                   self.locations[name].remove(loc2)
                   quit = 0
             elif loc2.startswith(location) or location.startswith(loc2):
                GVars.out.put('\nINFO: Parameter %s defined in both\n%s (%s default %s) and \n%s (%s default %s)'
                        %(name, loc2, typ2, val2, location, type, value),globals.INFO)
                GVars.out.put("Longer path wins, overriding the less specific instance.",globals.INFO)
                if location.startswith(loc2): # remove existing details
                   del self[(loc2,name)]
                   self.locations[name].remove(loc2)
                   quit = 0
             elif (
                (self[(loc2,name)]["RANGE"]==[] and range) or
                (self[(loc2,name)]["RANGE"]!=[] and not range)):
                GVars.out.put('\nINFO: Parameter %s defined in both\n%s (%s default %s) and \n%s (%s default %s)'
                        %(name, loc2, typ2, val2, location, type, value),globals.IMPINFO)
                GVars.out.put("Specification with range wins, overriding the less specific instance.",globals.IMPINFO)
                if range: # remove existing details
                   GVars.out.put("Value and range for %s %s set to %s [%s]."%(type,name,value,range),globals.IMPINFO)
                   del self[(loc2,name)]
                   self.locations[name].remove(loc2)
                   quit = 0
                else:
                   GVars.out.put("Value and range for %s %s set to %s %s."%(typ2,name,val2,self[(loc2,name)]["RANGE"]),globals.IMPINFO)
             else:
                GVars.out.put('\nWARNING: Parameter %s defined in both\n%s (default %s) and \n%s (default %s)'
                        %(name, loc2, val2, location, value),globals.WARN)
                GVars.out.put("Ignoring second instance",globals.WARN)
             GVars.out.pop()
             if quit: return
       if not self.warn:
          if self.has_key((location,name)): # same combo of both is there
             GVars.out.push()
             GVars.out.put('\nWARNING: Parameter %s is defined in %s multiple times. Ignoring new instance\n'
                           % (name,location),globals.WARN)
             GVars.out.pop()
             return
       if not self.typemap.has_key(type):
          raise SetupError("Invalid type %s for parameter %s" % (type,name))
       ans = { "TYPE":type, 
               "VALUE":value, 
               "CONST":const, 
               "RANGE": self.parseRange(range,type,name),
               "COMMENT":comment}
       self[(location,name)] = ans
       if not self.locations.has_key(name): self.locations[name] = []
       self.locations[name].append(location)

   # Takes a range object and returns a string describing the range object
   def printRange(self, rng):
       if not rng: return "Valid Values: Unconstrained"
       ans = []
       for a in rng:
           if type(a) == type("STRING"): 
              ans.append('"%s"' % a)
              continue
           min = a["min"]
           max = a["max"]
           if not min: min = "-INFTY"
           if not max: max = "INFTY"
           if min == max:
              ans.append("%s" % min)
           else:
              ans.append("%s to %s" % (min,max))
       return "Valid Values: %s" % ", ".join(ans)

   def writeRPInfo(self,fname=None,prefix=""):
       if not fname:
          fname = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.SetupParamsFilename)
       out = globals.IndentedOutput(4, open(fname, 'w'))
       rplist = [ (loc,name,a) for ((loc,name),a) in self.items() ]
       rplist.sort()
       currunit = None
       out.push()
       for (unitname,rpname,rpinfo) in rplist:
           if not unitname.startswith(prefix): continue
           # Starting a new unit
           if currunit != unitname:
              out.pop()
              out.put("")
              out.put(unitname)
              currunit = unitname
              out.push()
           if rpinfo["TYPE"] == "DOC":
              out.put(rpname)
           elif self.docIgnoreRe.match(rpinfo["COMMENT"]):
              continue
           else:
              if (rpinfo["CONST"]):
                 out.put('%s [%s] CONSTANT [%s]' % (rpname, rpinfo["TYPE"],rpinfo["VALUE"]))
              else:
                 out.put('%s [%s] [%s]' % (rpname, rpinfo["TYPE"],rpinfo["VALUE"]))
              if rpinfo["TYPE"] in ["INTEGER","REAL","STRING"]:
                 out.push()
                 out.put(self.printRange(rpinfo["RANGE"]))
                 out.pop()
           out.push()
           out.put(rpinfo["COMMENT"])
           out.pop()
       out.pop()
       out.file.close()
       
   def writeDuplications(self,fname):
       out = globals.IndentedOutput(4, open(fname,"w"))
       names = self.locations.keys()
       # for each name find the number of locations where its type is not DOC
       # if this number > 1, then consider it duplicated
       names = filter(lambda x: len(
                          filter(lambda y: self[(y,x)]["TYPE"] != "DOC", 
                                 self.locations[x])
                          ) > 1,names) # only multiply defined ones
       names.sort()
       for name in names:
           out.put(name + " defined in the following locations")
           out.push()
           for loc in self.locations[name]:
               if self[(loc,name)]["TYPE"] == "DOC": continue
               if self[(loc,name)]["CONST"]:
                  out.put("%s CONSTANT [%s]" % (loc,self[(loc,name)]["VALUE"]) )
               else:
                  out.put("%s [%s]" % (loc,self[(loc,name)]["VALUE"]) )
           out.pop()
           out.put("")
       # find names with no comments
       # names = all names with one definition
       badnames = []
       for (x,y) in self.locations.items():
           for loc in y:
               if not self[(loc,x)].get("COMMENT",None): badnames.append((x,loc))
       # badnames = all names where no COMMENT
       if badnames: 
          out.put("Runtime Parameters without any Comments")
          out.put("---------------------------------------")
       oldname = None
       for (name,loc) in badnames:
           if oldname != name:
              out.put("\n"+name)
              oldname = name
           out.push()
           out.put("no comment in %s/Config" % loc)
           out.pop()
       out.file.close()
                  
   # info about RuntimeParameter functions are hard coded into this function
   def genRulesCode(self,rpname,rptype,range):

       def clean(obj,typ):
           "add a . if obj consists only of numbers and we want it to be REAL"
           if typ == "INTEGER": return str(obj)
           rv = str(obj)
           if rv == "TINY": return "TINY(1.0)"
           if rv == "-TINY": return "-TINY(1.0)"
           suff = ""
           # check if obj ends with "+","-"
           # adjust the given numbers in that direction by a
           # small amount
           if rv[-1] in ["+","-"]:
              suff = rv[-1]
              rv = rv [:-1]
              suff = suff+"EPSILON(1.0)"
           if re.match("^[0-9]*$",rv): # only numbers found
              rv = rv + ".0"
           return rv+suff

       if not range: return ""
       if rptype not in ["STRING","INTEGER","REAL"]: return ""
       rv = ""; rvrest = ""
       if rptype == "STRING": # we are doing it for a string 
          numVals = len(range)
          mlen= 0
          for val in range:
              if type(val) != type("STRING"):
                 raise SetupError("'%s' is not a string value for Runtime Parameter '%s'" % (val,rpname))
              mlen = max(mlen,len(val))
          tpl = '"%%-%ds"' % mlen
          validValues = [ tpl % val for val in range ]
          rvrest = '  call rp_rules( "%s", %d, (/ %s /) )' % (rpname, numVals,", ".join(validValues))
          rv1 = ""
       elif rptype in ["INTEGER","REAL"]: # integers or reals
          # an example number of the given type
          if rptype == "INTEGER":
             typ= "1"
          else: typ= "1.0"
          numVals = len(range)
          minVals = []
          maxVals = []
          for rangespec in range:
              if rangespec["min"]:
                 minVals.append( clean(rangespec["min"],rptype))
              else: minVals.append( "-HUGE(%s)" % typ)
              if rangespec["max"]: 
                 maxVals.append( clean(rangespec["max"],rptype))
              else: maxVals.append( "HUGE(%s)" % typ)
          rvrest = '  call rp_rules( "%s", %d, (/ %s /), (/ %s /) )' % \
                   (rpname,numVals, ", ".join(minVals), ", ".join(maxVals))
       if (len(rvrest) > 0): 
          rv1 = ""
          while len(rvrest) > 132:
             rv1 = rv1 + rvrest[:130] + "&\n"
             rvrest = "&" + rvrest[130:]
          rv = rv1 + rvrest + "\n"
       return rv 

   # write out a default.par which has all default values commented
   def writeDefaultPar(self):
       header = ['Copy before editing!!!',
                 'Created by the setup script.',
                 '',
                 'Contains default values for all runtime parameters specific to this simulation', 
                 '']
       fname = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.RPDefaultParFilename)
       f = LazyFile(fname)
       f.write('## '+"\n## ".join(header)+"\n\n")
       locitems = self.locations.items()
       locitems.sort()
       for (rpname,rplocations) in locitems:
           if len(rplocations) != 1:
              GVars.out.put("%s found in multiple locations %s!! Programming Error?" % (rpname,rplocations))
           loc = rplocations[0]
           rpinfo = self[(loc,rpname)]
           # ignore documentation things and CONST runtime parameters (which cannot be changed by the user anyway)
           if rpinfo["TYPE"] in ["DOC","CONST"]: continue
           if rpinfo["TYPE"] == "BOOLEAN":
              rpval = ".%s." % rpinfo["VALUE"]
           else: rpval = rpinfo["VALUE"]
           f.write("# %-30s = %s\n" % (rpname,rpval))
       f.write("\n")
       f.close()


   # info about RuntimeParameter functions are hard coded into this function
   def writeCode(self, configInfo):

       def getExpo(value):
          # extracts a conservative estimate of the exponent from a real constant,
          # this is not perfect but will do for the actual purpose; estimate whether
          # a numerical runtime parameter value from a Config file is so large that
          # it may be outside the range of exponents for 4-byte reals.
           if not value: return 0
           m = self.numExpoRe.match(value)
           if m:
              expo = m.group("expo")
              try:
                 before = len(m.group("before"))
              except:
                 before = 0
              try:
                 before2 = len(m.group("before2"))
              except:
                 before2 = 0
              b=max(before,before2)
              a = 0
              if b==0:
                 try:
                    a = len(m.group("aft0"))
                 except:
                    a = 0
              try:
                 exp = int(expo)
                 ex = max(exp, exp + b - 1)
                 ex0 = min(exp, exp + b - 1 - a)
                 return max(abs(ex),abs(ex0))
              except:
                 return 0
           else: return 0

       header = '! Runtime-settable parameter initializations;\n'\
                '! generated by setup script.\n'\
                '! Do not edit!\n\n'\
                'subroutine rp_initParameters(parmfile)\n\n'\
                'character(len=*) :: parmfile\n'\
                'call rp_createParmList()\n'\
                'return\n'\
                'end subroutine rp_initParameters\n\n\n'\
                'subroutine rp_createParmList ()\n\n'\
                'use RuntimeParameters_interface, ONLY : RuntimeParameters_add\n'\
                'use RuntimeParameters_data, ONLY : TYPE_CONST, TYPE_VAR\n\n'\
                'implicit none\n\n'\
                'integer,parameter :: r8=kind(1.0)\n'\
                '#include "constants.h"\n'\
                '#include "rp_rules.h"\n\n'

       fname = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.RPInitParamsFilename)
       f = LazyFile(fname)

       f.write(header)

       locationslist = self.locations.items()
       locationslist.sort()
       
       for (rpname,rplocations) in locationslist:
           if len(rplocations) != 1:
              GVars.out.put("%s found in multiple locations %s!! Programming Error?" % (rpname,rplocations))
           loc = rplocations[0]
           rpinfo = self[(loc,rpname)]
           # Ignore "parameters" which are documentation
           if rpinfo["TYPE"] == "DOC": continue
           if rpinfo["TYPE"] == "BOOLEAN": 
              rpvalue = ".%s."
           else: rpvalue = "%s"


           
           rpvalue = rpvalue % rpinfo["VALUE"]

           #print rpvalue
           if rpinfo["TYPE"] == "REAL":
             if rpvalue == "TINY":
               rpvalue = "TINY(1.0)"
             elif rpvalue == "HUGE":
               rpvalue = "HUGE(1.0)"
             elif rpvalue == "-HUGE":
               rpvalue = "-HUGE(1.0)"
             elif getExpo(rpvalue) > 37:
               rpvalue = rpvalue.upper().replace('D','e') + "_r8"
           #elif rpinfo["TYPE"] == "INTEGER":
            # if rpvalue == "MAXINT"
             #   rpvalue = "
           
           if rpinfo["CONST"]:
              f.write('  call RuntimeParameters_add( "%s", %s, %s)\n'%(rpname, rpvalue, "TYPE_CONST"))
           else:
              f.write('  call RuntimeParameters_add( "%s", %s)\n'%(rpname, rpvalue))
           # The code for interval and enum checking should come here
           rules = self.genRulesCode(rpname,rpinfo["TYPE"],rpinfo["RANGE"])
           if rules: f.write(rules)
           
#        #Allocate enough io_plot_var names for all variables in the simulation.
#        #Not needed here any more, now instead calling addRP from getRPInfo (in unitUtils.py) for these.
#        numRegularPlotVars = configInfo['max_plot_vars']
#        for i in range(1,numRegularPlotVars+1):
#            f.write('  call RuntimeParameters_add( "%s", %s)\n'%(('plot_var_' + '%d'%i), '"none"'))

       f.write('\nreturn\nend subroutine rp_createParmList\n\n')
       f.close()

