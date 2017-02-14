"""Code related to Variable COMMENTS"""
import string
import re
import sys
import UserDict
import os.path

# relative imports needed
from . import setup_globals
from .setup_globals import gvars, SetupError


class VarInfo(UserDict.UserDict):
    """Class which keeps track of variables and related info
    To each variable we keep track of the following info
   
    "NAME" of the variable
    "ATTRIBUTES" - Conservation laws satisfied
    "COMMENT" - Info regarding use of variable
   
    self is a dictionary mapping ("LOCATION","NAME","TYPE") to other info
    """

    def addVar(self,name,type="Unknown",attribs=[],location="",comment=""):
       if self.has_key((location,name,type)): # same combo of both is there
          gvars.out.push()
          gvars.out.put('\nWARNING: %s %s is defined in %s multiple times. Ignoring new instance\n'
                        % (type,name,location),setup_globals.WARN)
          gvars.out.pop()
          return
       ans = { "ATTRIBS":attribs, "COMMENT":comment}
       self[(location,name,type)] = ans

    def write_var_info(self,fname=None,prefix=""):
       if not fname:
          #fname = os.path.join(gvars.flash_src_dir,gvars.project_build_dir, setup_globals.SETUP_VARS_FILENAME)
          fname = os.path.join(gvars.project_setup_dir, setup_globals.SETUP_VARS_FILENAME)
       out = setup_globals.IndentedOutput(4, open(fname, 'w'))
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
       
