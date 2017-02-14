#!/usr/bin/env python

"""This program parses all the .c files in the current directory and finds the header files
each file depends on. Having done that, it appends those dependencies to Makefile.Depends."""

import re, os, os.path, string, sys

incpatt = """^\s*[#]\s*include\s+["']?([A-Za-z._0-9]+)['"]?\s+.*"""

# searches for file "name" in a specified list of directories 
# and returns the absolute name
def findINCmatch(name):
    global INCDirs
    global INCcache
    if INCcache.has_key(name): return INCcache[name]
    for d in INCDirs:
        dname = os.path.join(d,name)
        if os.path.isfile(dname):
           INCcache[name] = dname
           return dname
    INCcache[name]=""
    return ""
 
# Given a filename computes the include files it depends on
# returns the include files as a list of strings
def depends(filename):
    """Handle one file"""
    incs = {}
    if not os.path.isfile(filename): return []
    for x in file(filename).readlines():
        m = incRE.match(x)
        if m: 
           incs[m.group(1)] = 1
           continue
    return [findINCmatch(x) for x in incs.keys()]


##################
def main():
   global incRE
   global INCcache
   global filelist, INCDirs

   # setup global variables
   incRE = re.compile(incpatt,re.I)

   # Process all the files 
   INCcache = {}

   filelist = os.listdir(".") # All files in current directory
   INCDirs = [""] # List of directories which contain include files. "" is current directory
   toProcess = [] # List of files to generate dependency info for

   # Process argument list and isolate filenames , inc dirs, arbitrary flags 
   for arg in sys.argv[1:]: # for each option given
       if not arg: continue
       if arg[0] == "-" and len(arg) > 2 and arg[1] == "I": # option -I
          INCDirs.append(arg[2:])
       elif arg[0] != "-": # not an option so must be a filename to process
          if arg[0] != '*': toProcess.append(arg) # unprocessed wild card means nu such files
       else: 
          pass # Ignore all other options

   ofd = file("Makefile.Depend","a")
   ofd.write("\n###Automatically appended dependencies of C files on headers###\n")
   for filename in toProcess:
       basename,ext = os.path.splitext(filename)
       incs = depends(filename)
       for h in incs: # Add more files to process if required
           if (h in filelist) and (h not in toProcess): # new header file in current directory
              toProcess.append(h)
       incpart  = " ".join(incs)
       ofd.write("%s.o: %s\n" % (basename, incpart))

   ofd.close()
   

if __name__=="__main__":
   main()

