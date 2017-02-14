
# code related to Variable COMMENTS and stuff

__all__ = [ "VarInfo", 
          ]

##############################################

import globals
from globals import *   # GVars and SetupError

import string, re, sys, UserDict, os.path

###############################################################

class VarInfo(UserDict.UserDict):
   # Class which keeps track of variables and related info
   # To each variable we keep track of the following info
   #
   # "NAME" of the variable
   # "ATTRIBUTES" - Conservation laws satisfied
   # "COMMENT" - Info regarding use of variable
   # 
   # self is a dictionary mapping ("LOCATION","NAME","TYPE") to other info

   def addVar(self,name,type="Unknown",attribs=[],location="",comment=""):
       if self.has_key((location,name,type)): # same combo of both is there
          GVars.out.push()
          GVars.out.put('\nWARNING: %s %s is defined in %s multiple times. Ignoring new instance\n'
                        % (type,name,location),globals.WARN)
          GVars.out.pop()
          return
       ans = { "ATTRIBS":attribs, "COMMENT":comment}
       self[(location,name,type)] = ans

   def writeVarInfo(self,fname=None,prefix=""):
       if not fname:
          fname = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.SetupVarsFilename)
       out = globals.IndentedOutput(4, open(fname, 'w'))
       varlist = [ (name,type,loc,a) for ((loc,name,type),a) in self.items() ]
       varlist.sort()
       currunit = None
       for (varname,vartype,unitname,varinfo) in varlist:
           if not unitname.startswith(prefix): continue
           out.put('Name: %s' % varname)
           out.put('Location: %s' % unitname)
           out.put('Type: %s' % vartype)
           out.put('Description:')
           out.push()
           out.put(varinfo["COMMENT"])
           out.pop()
           out.put("")
       out.file.close()
       
