
# code to generate files, based on a template

__all__ = [ "Template",
          ]

##############################################

from globals import *
import globals
from lazyFile import * # for LazyFile class

import string, re, os.path, UserDict

###############################################################

class Template(UserDict.UserDict):
# generates file filename using template as the template
# Lines starting with comment in the template are suppressed
# other lines are % substituted against the given
# dictionary
# returns if generatedFile is same as existing file

  def __init__(self,template,comment="##"):
      UserDict.UserDict.__init__(self)
      if not os.path.isfile(template):
         raise SetupError("Template file %s not found" % template)
      self.template = os.path.abspath(template)
      self.comment = comment
      self.check = None # Dont perform checks on what is being accessed
      self.__getitem = UserDict.UserDict.__getitem__
      self.__setitem = UserDict.UserDict.__setitem__

  def __setitem__(self,key,value):
      UserDict.UserDict.__setitem__(self,key,value)
      if type(value) in [type([]),type(()),type({})]:
         self.__setitem(self,"COUNT_"+key,len(value))
         try:
           maxlen = max([len(str(x)) for x in value])
         except ValueError:
           maxlen = 0
         self.__setitem(self,"MAXLEN_"+key, maxlen)

  def update(self,dict):
      for (k,v) in dict.items(): self[k] = v

  # Fancy printing of list entries
  def __getitem__(self,key):
      parts = key.split("|")
      if len(parts) == 1: # No | found
         if key.find("!") < 0: #simple variables
            self.__checksimple(key)
            return self.__getitem(self,key)
         else: return self.__printlist(key,"")
      # we have a nil value
      var = parts[0] # variable to print
      nil = parts[1] # print if NULL or empty or 0
      if var.find("!") < 0: # we are printing a string or integer
         self.__checksimple(var)
         if self[var]:
            return self.__getitem(self,var)
         else: return nil
      # Now we need to handle sequence type with a nil value
      return self.__printlist(var,nil)

  # Check if value is a simple datatype
  def __checksimple(self,key):
      if not self.check: return 
      if not self.has_key(key):
         raise SetupError("Variable %s unknown in template %s" 
                    % (key,self.template),globals.IMPINFO)
      tkey = type(self.__getitem(self,key))
      if tkey in [type(0),type(""),type(None)]: return
      GVars.out.put("WARNING: In template %s you are using the variable '%s' as a simple variable" 
                    % (self.template,key),globals.IMPINFO)
      GVars.out.push()
      GVars.out.put("But it has type %s" % tkey,globals.IMPINFO)
      GVars.out.pop()

  # Allow fancy stuff inside %(...)s
  # Never use self[...] in this code as we are implementing the code for self[...]
  def __printlist(self,key,nilvalue):
      parts = key.split("!")
      if len(parts) > 3:
         raise SetupError("Bad Template Syntax (%s) in %s\n Found more than 2 !"% (key,self.template))
      info = { "num" : 1,  # block size 
               "var" : "", # name of list variable
               "in"  : "", # inner separator
               "out" : "", # outer separator
               "rpad": "", # pad all entries with spaces to right so they come out as same length
             }
      mobj = re.match("^(?P<num>[0-9]*)(?P<var>[a-zA-Z_]*)(?P<rpad> *)$",
                      parts[0])
      if not mobj:
         raise SetupError("Bad template syntax (%s) in %s" % (key,self.template))
      info.update(mobj.groupdict())
      if not self.has_key(info["var"]):
         raise SetupError("Unknown variable %s while expanding template %s" % (info["var"],self.templte))
      vals = self.__getitem(self,info["var"])
      if type(vals) not in [type([]),type(())]:
         raise SetupError("Variable %s has type %s in template %s -- not a sequence type" 
                          % (info["var"],type(vals),self.template))
      # Handle trivial case - empty sequence
      if not vals:
         return nilvalue
      # compute block size
      if info["num"]:
         info["num"] = int(info["num"])
      else: info["num"] = 1
      # parse meaning of parts
      if info["num"] == 1:
         # pretend block size 2 with inner and outer sep same
         info["num"] = 2
         info["in"] = parts[1]
         info["out"] = parts[1]
      else:
         if len(parts) == 2:
            raise SetupError("""Bad Template Syntax (%s) in %s
    For non-trivial block size you must specify inner and outer separators""" % (key,self.template))
         info["in"] = parts[1]
         info["out"]= parts[2]
      # form the output string
      count = 0
      ans = []
      if info["rpad"]:
         tstr = "%%-%ds"% self.__getitem(self,"MAXLEN_"+info["var"])
      else: tstr = "%s"
      for val in vals:
          count += 1
          ans.append(tstr % val)
          if count % info["num"] == 0: # end of block
             ans.append(info["out"])
          else: ans.append(info["in"])
      # we have an extra separator at the end
      return string.join(ans[:-1],"")
      
  def generate(self,filename):

    def repl(mobj):
        s = mobj.group(0)
        if s == r"\\": return '\\'
        elif s == r"\t": return '\t'
        elif s == r"\r": return '\r'
        elif s == r"\n": return '\n'
        else: return s
 
    infd = open(self.template,"r")
    # if the out file already exists and is a link, kill it
    if os.path.islink(filename): os.remove(filename)
    outfd= LazyFile(filename)
    length = len(self.comment)
    self.check = 1 # Perform checks on what is being accessed
    for line in infd.readlines():
        if line[:length] == self.comment: continue # skip comments
        # perform % substitution and replace \{\rtn} sequences
        outfd.write( re.sub(r"\\.",repl,line) % self)
    infd.close()
    outfd.close()
    self.check = None # disable checks
    return outfd.samefile

