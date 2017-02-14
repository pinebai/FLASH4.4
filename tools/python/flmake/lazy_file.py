"""Implements a Lazy File class

Associated with a file which is regenerated afresh and never read from.
This class doesn't write the actual file till the file is closed
even then it writes only if the new contents and the current contents 
are different.
"""

import os.path

# relative imports needed!    
from . import setup_globals
from .setup_globals import gvars, SetupError


class LazyFile:
   """Defines a file like class which supports the following:

    * __init__ associates a real file with this class
    * write adds a line for writing to the file (but does not actually write)
    * close reads the real file, compares it with buffer contents
      if doesnot match rewrites the entire file, else silently returns.

    NOTE: Useful only if you want to rewrite the entire file, small files (memory requirement)

    Usage:: 

        file = LazyFile("filename")
        file.write(string1)
        file.write(string2)
        ...
        file.write(stringn)
        file.close() --> if reqd this is the one which really does all the writing

   """

   def __init__(self,filename):
       self.filename = os.path.abspath(filename)
       self.contents = []
       self.curr = "" # the current partial line

   def write(self,what):
       if not what: return # Nothing to write
       a = what.split("\n") # if multiple lines are in what
       if self.curr:
          self.contents.append("%s%s\n" % (self.curr,a[0]))
          self.curr = ""
          del a[0]
       self.curr = a[-1]
       del a[-1]
       self.contents.extend(["%s\n"%x for x in a])

   # Return 1 if new contents match old contents else return None
   def match(self):
       try:
         myfile = open(self.filename,"r")
       except:
         return None
       ctr = 0
       try:
         for x in myfile.xreadlines():
           if x != self.contents[ctr]:
              myflie.close()
              return None
           else: ctr += 1
         # writing more lines than before
         if len(self.contents[ctr:]) > 0: 
            myfile.close()
            return None
       except:
         # new contente is shorter than old contents
         myfile.close()  
         return None
       myfile.close()
       return 1
           
   def close(self):
       # is there is a partial line
       if self.curr: 
          self.contents.append(self.curr)
          self.curr = ""
       # no need to rewrite the file
       self.samefile = self.match()
       if self.samefile==1: 
          gvars.out.put("file %s did not really change" % self.filename, setup_globals.INFO)
          return
       myfile = open(self.filename, "w")
       for x in self.contents: myfile.write(x)
       myfile.close()

