

# implementation of LinkFileList class

__all__ = ["LinkFileList"]

__doc__ = """
This class contains the list of files to be linked to the object directory
 
 list = LinkFileList(objectDir)
 list.cleanObjectDir()
   ---- after all units have been generated ----
 list.getDontLinkList(units,linkifUDPairs)
 list.linkFiles (sourcedir1)
 list.linkFiles (sourcedir2)
 ......
 list.linkFiles (sourcedir9)
 list.reallyLink() -- this is the one which really makes the links
"""

#################################

import globals
from globals import * # GVars and SetupError
from utils import * # for GetRelPath

import types, os, string, shutil
try:
  import zlib
except ImportError, e:
  # print ("Tried importing zlib, %s, continuing anyway..." % e)
  #  GVars.out.put doesn't seem to work here - KW
  pass

class LinkFileList:
  
  # Get all user options, Global Variables and all units
  def __init__(self,objDir):
    self.objDir = objDir
    # In case of noClobber dont touch files with following extension initially
    self.noClobberExt = ['.so', '.a', '.o', '.unit', '.F90', '.c', '.fh', 
                    ".F",".h",".mod",'.py'] 
    # extensions of files are never deleted (unless the source is deleted)
    self.derExt = [".o",".mod"]
    # precedence for extensions in each group
    groups = [['.f90', '.F90', '.f', '.c', '.F', '.C', '.cxx'], '.py', '.fh', ['.h','.hxx']]
    #put extensions that aren't grouped into a one-item list (save typing)
    self.groupdict = {} # dictionary mapping an extension to list of related extensions
    self.exts = []
    self.groups = []
    # DEV as far as I can tell, neither 'self.groups' nor
    # 'self.groupdict' (see below) is ever used for anything
    # in this file or in any file that instantiates this class.
    #  - nttaylor
    for x in groups:
        if type(x) == types.ListType:
           self.groups.append(x)
           self.exts.extend(x)
           for y in x: self.groupdict[y] = x
        else: 
           self.groups.append([x])
           self.exts.append(x)
           self.groupdict[x] = [x]

    self.links = {} # dictionary mapping linked name to real name (absolute path)

  def addLink(self,realname,linkname):
     # By storing this info in a dictionary we have that later addLinks
     # overwrite earlier addlinks with the same name
     # checkLinks method checks for abc.c and abc.F90 type issues
     self.links[linkname] = realname

  def isUnitPresent(self,uname,unitList):
     """Check if uname or one of its children is present in units"""
     for unitname in unitList:
       if unitname.startswith(uname): return 1
     return None


  def doLINKIFOverrides(self, unitlist, udpairs):
    # unitlist = list of names of all units used in this simulation
    # udpairs = list of all (fname,uname) pairs found in "LINKIF" statements.
    functionOverrideList = []
    allLinkifFunctions = {}

    # make dictionary mapping "LINKIF" function names
    # to lists of all units on which that function depends.
    for (fname, uname) in udpairs:
      if allLinkifFunctions.has_key(fname):
        allLinkifFunctions[fname].append(uname)
      else:
        allLinkifFunctions[fname] = [uname]

    # If each and every unit on which a function depends is
    # included, add that function to 'functionOverrideList'
    for fname in allLinkifFunctions.keys():
      for uname in allLinkifFunctions[fname]:
        if not self.isUnitPresent(uname, unitlist):
          # One of the units on which this function depends is not
          # present, so break and skip the append statement below.
          break
      else:
        functionOverrideList.append(fname)

    # Do the overrides last.  The later call to 'addLink'
    # ensures that these values will override any previous ones
    functionOverrideList.sort()
    functionOverrideList.reverse()
    for file in functionOverrideList:
        # Pull out the parts delimited by "."
        # throw out all but the first 2, and join them back up
        # file = "a/b/c/d.e.f.g" -> outfile = "d.e"
        outfile = string.join( string.split(os.path.basename(file),".")[:2],".")
        self.addLink(os.path.join(GVars.sourceDir, file), outfile)


  # Input: unitlist = list of names of all units
  #        udpairs = list of all (fname,uname) pairs found
  #  Compute list of files which should not be linked to
  #  based on the LINKIF directives
  # (fname,uname) in LINKIF means
  #   link in fname only if uname is used as a UNIT, otherwise dont link in fname
  # Algo
  # ans = all fnames occurring in LINKIF
  # if (fname,uname) in LINKIF and unit in units such that uname is a 
  # prefix of unit.name, remove fname from ans
  def getDontLinkList(self,unitlist,udpairs):
     """Return a list of files which should not be linked to"""
     fdict = {}
     deldict = {} # List of file names which should not be linked
     for (fname,uname) in udpairs: 
        deldict[fname] = 1 
        try:
           fdict[fname].append(uname)
        except KeyError: 
           fdict[fname] = [uname]
     # now fdict is a dictionary mapping filenames to units they depend on
     # If a file depends on MULTIPLE units, it will be included only if ALL
     # the units it depends on are in the units list
     for fname in deldict.keys():
        # find num units in fdict[fname] not in unitlist
        ans = [x for x in fdict[fname] if not self.isUnitPresent(x,unitlist)]
        if not ans: del deldict[fname]
     # Return the absolute path of file names left in deldict
     self.DontLink = [os.path.join(GVars.sourceDir,x) for x in deldict.keys()]
     GVars.out.put("%d Files which will not be linked:"%len(self.DontLink),globals.DEBUG)
     GVars.out.push()
     for m in self.DontLink: GVars.out.put(m,globals.DEBUG)
     GVars.out.pop()


  def linkFiles(self, fromdir):
    """
    Link files in with right extension in fromdir to object directory

    Note: if blah.c is in fromdir and blah.F90 is in the current directory,
    blah.F90 gets deleted. Same thing for other groups of extensions.

    We now also accept files like abc.F90.x.y.z to stand for .F90 files. When 
    this file is linked it will be linked as "abc.F90" The real extension of a file
    is the string between the first and the second dots, (or string after the first,
    if there is no second dot)

    NOTE: We dont actually link the files, only queue them up for linking. The real
    linking is done by calling reallyLink method
    """
    
    files = []
    GVars.out.push()
    for file in os.listdir(fromdir):
        # DEV come back here to re-insert fix - ntt
        parts = os.path.basename(file).split(".",1)
        if len(parts) > 1:
           ext = "." + parts[1]
        else: ext = None
        if ext in self.exts:
           fullname = os.path.join(fromdir,file)
           if fullname in self.DontLink:
              GVars.out.put("Not linking %s" % fullname)
           else: files.append(fullname)
    GVars.out.pop()

    #sort to get consistent behavior when 2 files with equivalent extensions
    #are in fromdir    
    files.sort()
    files.reverse() 
    for file in files:
        # Pull out the parts delimited by "."
        # throw out all but the first 2, and join them back up
        # file = "a/b/c/d.e.f.g" -> outfile = "d.e"
        outfile = string.join( string.split(os.path.basename(file),".")[:2],".")
        self.addLink(file, outfile)


  def reallyLink(self):
    """Function serves to link files whether portable is selected or not"""
    GVars.out.put("Really linking files into object dir",globals.INFO)
    # ensure we are not going to link abc.c and abc.F90 in the current run
    self.checkLinks() 
    if GVars.noClobber: self.changedLinks()
    pwd = os.getcwd()
    # run through the links and do the appropriate thing
    if GVars.portable:
       for (link,real) in self.links.items():
          if os.path.isfile(link): os.remove(link)
          shutil.copy2(real,link)
    else:
       for (link,real) in self.links.items():
          if os.path.isfile(link): os.remove(link)
          else:
#            print("Not a file: %s\n" % link)
            if not os.path.exists(link) and  os.path.islink(link) and os.path.lexists(link):
              GVars.out.put("Seems to be a broken symlink, will remove: %s" % link,globals.DEBUG)
              os.remove(link)
          try:
            os.symlink(getRelPath(real,pwd),link)
          except OSError:
            try:
              relpath = getRelPath(real,pwd)
            except:
              GVars.out.put("Failure to get path of %s relative to %s for symlinking\n" % (real,pwd),globals.IMPINFO)
              raise
            GVars.out.put("Failure symlinking %s to %s\n" % (relpath,real),globals.IMPINFO)
            raise

  # since this can be called in a portable setting, we also need to handle the case
  # that there is a abc.F from the previous run and an abc.F90 from this run
  # in this case, we need to remove the abc.F from the prev run as well
  def rmLink(self,linkname):
    """Remove the file linkname and associated files"""
    GVars.out.put("Removing %s and associates" % linkname, globals.DEBUG);
    linkbase,ext = os.path.splitext(linkname)
    for exta in self.groupdict.get(ext,[ext]): # all related extensions
        try: os.remove(linkbase+exta)
        except: pass #if file does not exist ignore the error
    try: os.remove(linkbase+".o")
    except OSError: pass
    try: os.remove(linkbase.lower()+".mod")
    except OSError:
      try: os.remove(linkbase.upper()+".mod")
      except OSError: pass

  # return if fname is in noClobberExceptionList
  def isSetupFile(self,fname):
      for prefix in globals.noClobberExceptionList:
          if fname[:len(prefix)] == prefix: return 1
      return 0 

  def cleanObjectDir(self):

      if not os.path.isdir(self.objDir): os.makedirs(self.objDir)

      if GVars.noClobber: #dont delete all files
         GVars.out.put('removing files except source and object files (if any)',globals.INFO)
      else:
         GVars.out.put('removing old links in build directory %s'%self.objDir,globals.INFO)

      for fname in os.listdir(self.objDir):
          # ignore directories
          if os.path.isdir(os.path.join(self.objDir, fname)): continue
          # what extension?
          if (GVars.noClobber and (os.path.splitext(fname)[1] in self.noClobberExt)): 
             continue
          else:
             # should we keep this file
             if not self.isSetupFile(fname) or (os.path.splitext(fname)[1] == '.o'):
                os.remove(os.path.join(self.objDir, fname))

  def checkLinks(self):
    # check to see if we are trying to link in a ABC.F90 as well as an ABC.C, if so 
    # inform user of potential problem and quit
    keys = self.links.keys()
    for file in keys:
        base,ext = os.path.splitext(file)
        be = base+ext
        for ext2 in self.groupdict[ext]:
            be2 = base+ext2
            if ext2 != ext and be2 in keys:
               raise SetupError('ERROR Checking Links: Both %s and %s should not be in object directory' %
                                (be,be2))

  # The following method should check if the two files are the same by computing the adler32 checksum.
  # (It is theoretically possibile but UNLIKELY that two different text files produce the same checksum.)
  # However, this version of linkFiles.py is adapted for machines where the zlib module
  # of Python is unavailable, so if we cannot use it (as indicated by a NameError exception),
  # let the files appear different.
  def samefile(self,filea,fileb):
      try:
        fa = open(filea)
        cksuma = zlib.adler32(fa.read())
        fa.close()
      except IOError:
        cksuma = None
      except NameError:
        fa.close()
        cksuma = "foo"

      # non-existant file is not equal to any file
      if not cksuma: return None

      try:
        fb = open(fileb)
        cksumb = zlib.adler32(fb.read())
        fb.close()
      except IOError:
        cksumb = None
      except NameError:
        fb.close()
        cksumb = "bar"

      return cksuma == cksumb

  def changedLinks(self):
    # Find the list of files in the current dir with significant extension
    # compare with self.links
    # links no longer present and are not derived files get deleted
    # for links where target has changed the .o file gets deleted
    # target same and has been updated is no problem (make takes care of it)
    for fname in os.listdir("."):
        # Dont touch files which should be retained between runs
        # or handled elsewhere in the code
        if self.isSetupFile(fname): continue
        # Do not touch derived files  (.o,.mod)
        if os.path.splitext(fname)[1] in self.derExt: continue
        # the place from where the file should be copied (for the new run)
        # this is None if file not required for new run
        newtarget = self.links.get(fname,None)
        if newtarget:
           if os.path.islink(fname): # fname is a symlink, so check if it points to new target
              same = ( os.path.abspath(os.readlink(fname)) == newtarget)
           else: 
              same = self.samefile(fname, newtarget)
        else: 
              same = None
        # target changed or the link is no longer required
        if not same:
           GVars.out.put("current file [%s] is different from [%s]. Removing it"% (fname,newtarget),globals.INFO)
           self.rmLink(fname)

