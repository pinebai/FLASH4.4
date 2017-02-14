"""Implementation of LinkFileList class."""

import os
import shutil

try:
    import zlib
except ImportError as e:
    # print ("Tried importing zlib, %s, continuing anyway..." % e)
    # gvars.out.put doesn't seem to work here - KW
    pass

from . import setup_globals
from .setup_globals import gvars, SetupError
from .utils import get_rel_path


class LinkFileList(object):
    """This class contains the list of files to be linked to the object directory

    list = LinkFileList(project_setup_dir)
    list.clean_setup_dir()
    ---- after all units have been generated ----
    list.get_dont_link_list(units,linkifUDPairs)
    list.link_files (sourcedir1)
    list.link_files (sourcedir2)
    ......
    list.link_files (sourcedir9)
    list.really_link() -- this is the one which really makes the links
    """

    def __init__(self, setup_dir):
        """Get all user options, Global Variables and all units."""
        self.setup_dir = setup_dir

        # In case of noClobber dont touch files with following extension initially
        self.no_clobber_ext = ['.so', '.a', '.o', '.unit', '.F90', '.c', '.fh',
                               ".F", ".h", ".mod", '.py']

        # extensions of files are never deleted (unless the source is deleted)
        self.der_ext = [".o", ".mod"]

        # precedence for extensions in each group
        groups = [['.f90', '.F90', '.f', '.c', '.F', '.C'], '.py', '.fh', '.h']

        # put extensions that aren't grouped into a one-item list (save typing)
        # dictionary mapping an extension to list of related extensions
        self.groupdict = {}
        self.exts = []
        self.groups = []

        # DEV as far as I can tell, neither 'self.groups' nor
        # 'self.groupdict' (see below) is ever used for anything
        # in this file or in any file that instantiates this class.
        #  - nttaylor
        for x in groups:
            if isinstance(x, list):
                self.groups.append(x)
                self.exts.extend(x)
                for y in x:
                    self.groupdict[y] = x
            else:
                self.groups.append([x])
                self.exts.append(x)
                self.groupdict[x] = [x]

        # dictionary mapping linked name to real name (absolute path)
        self.links = {}

    def add_link(self, realname, linkname):
        """By storing this info in a dictionary we have that later add_links
        overwrite earlier addlinks with the same name
        check_links method checks for abc.c and abc.F90 type issues"""
        self.links[linkname] = realname

    def is_unit_present(self, uname, unit_list):
        """Check if uname or one of its children is present in units"""
        for unitname in unit_list:
            if unitname.startswith(uname):
                return 1
        return None

    def do_linkif_overrides(self, unitlist, udpairs):
        # unitlist = list of names of all units used in this simulation
        # udpairs = list of all (fname,uname) pairs found in "LINKIF" statements.
        function_override_list = []
        all_linkif_functions = {}

        # make dictionary mapping "LINKIF" function names
        # to lists of all units on which that function depends.
        for (fname, uname) in udpairs:
            if fname in all_linkif_functions:
                all_linkif_functions[fname].append(uname)
            else:
                all_linkif_functions[fname] = [uname]

        # If each and every unit on which a function depends is
        # included, add that function to 'function_override_list'
        for fname in all_linkif_functions.keys():
            for uname in all_linkif_functions[fname]:
                if not self.is_unit_present(uname, unitlist):
                    # One of the units on which this function depends is not
                    # present, so break and skip the append statement below.
                    break
            else:
                function_override_list.append(fname)

        # Do the overrides last.  The later call to 'add_link'
        # ensures that these values will override any previous ones
        function_override_list.sort()
        function_override_list.reverse()
        for file in function_override_list:
            # Pull out the parts delimited by "."
            # throw out all but the first 2, and join them back up
            # file = "a/b/c/d.e.f.g" -> outfile = "d.e"
            outfile = ".".join(os.path.basename(file).split(".")[:2])
            self.add_link(os.path.join(gvars.source_dir, file), 
                          os.path.join(gvars.project_setup_dir, outfile))

    def get_dont_link_list(self, unitlist, udpairs):
        """Return a list of files which should not be linked to

        Input: unitlist = list of names of all units
        udpairs = list of all (fname,uname) pairs found
        Compute list of files which should not be linked to
        based on the LINKIF directives
        (fname,uname) in LINKIF means
        link in fname only if uname is used as a UNIT, otherwise dont link in fname
        Algo
        ans = all fnames occurring in LINKIF
        if (fname, uname) in LINKIF and unit in units such that uname is a
        prefix of unit.name, remove fname from ans
        """
        fdict = {}
        deldict = {}

        # List of file names which should not be linked
        for (fname, uname) in udpairs:
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
            ans = [x for x in fdict[fname] if not self.is_unit_present(x, unitlist)]
            if not ans:
                del deldict[fname]

        # Return the absolute path of file names left in deldict
        self.dont_link = [os.path.join(gvars.source_dir, x) for x in deldict.keys()]
        gvars.out.put("{0} Files which will not be linked:".format(len(self.dont_link)),
                      setup_globals.DEBUG)
        gvars.out.push()
        for m in self.dont_link:
            gvars.out.put(m, setup_globals.DEBUG)
        gvars.out.pop()

    def link_files(self, fromdir):
        """Link files in with right extension in fromdir to object directory

        Note: if blah.c is in fromdir and blah.F90 is in the current directory,
        blah.F90 gets deleted. Same thing for other groups of extensions.

        We now also accept files like abc.F90.x.y.z to stand for .F90 files. When
        this file is linked it will be linked as "abc.F90" The real extension of a file
        is the string between the first and the second dots, (or string after the first,
        if there is no second dot)

        NOTE: We dont actually link the files, only queue them up for linking. The real
        linking is done by calling really_link method
        """
        files = []
        gvars.out.push()
        for file in os.listdir(fromdir):
            # DEV come back here to re-insert fix - ntt
            parts = os.path.basename(file).split(".", 1)
            if len(parts) > 1:
                ext = "." + parts[1]
            else:
                ext = None

            if ext in self.exts:
                fullname = os.path.join(fromdir, file)
                if fullname in self.dont_link:
                    gvars.out.put("Not linking {0}".format(fullname))
                else:
                    files.append(fullname)
        gvars.out.pop()

        # sort to get consistent behavior when 2 files with equivalent extensions
        # are in fromdir
        files.sort()
        files.reverse()
        for file in files:
            # Pull out the parts delimited by "."
            # throw out all but the first 2, and join them back up
            # file = "a/b/c/d.e.f.g" -> outfile = "d.e"
            outfile = ".".join(os.path.basename(file).split(".")[:2])
            outfile = os.path.join(gvars.project_setup_dir, outfile)
            self.add_link(file, outfile)

    def really_link(self):
        """Function serves to link files whether portable is selected or not"""
        gvars.out.put("Really linking files into object dir", setup_globals.INFO)

        # ensure we are not going to link abc.c and abc.F90 in the current run
        self.check_links()
        if gvars.noClobber:
            self.changed_links()
        pwd = os.getcwd()

        # run through the links and do the appropriate thing
        if gvars.portable:
            for (link, real) in self.links.items():
                if os.path.isfile(link):
                    os.remove(link)
                shutil.copy2(real, link)
        else:
            for (link, real) in self.links.items():
                if os.path.isfile(link):
                    os.remove(link)
                else:
                    #print("Not a file: %s\n" % link)
                    if not os.path.exists(link) and  os.path.islink(link) \
                       and os.path.lexists(link):
                        gvars.out.put("Seems to be a broken symlink, will remove: {0}".format(link),
                                      setup_globals.DEBUG)
                        os.remove(link)
                try:
                    # For absolute symlinks
                    os.symlink(os.path.abspath(real), link)
                    # For relative symlinks
                    #os.symlink(get_rel_path(real,pwd),link)
                except OSError:
                    gvars.out.put("Failure symlinking {0} to {1}\n".format(get_rel_path(real, pwd),
                                  real) % link, setup_globals.IMPINFO)
                    raise

    # since this can be called in a portable setting, we also need to handle the case
    # that there is a abc.F from the previous run and an abc.F90 from this run
    # in this case, we need to remove the abc.F from the prev run as well
    def rm_link(self, linkname):
        """Remove the file linkname and associated files"""
        #print "Removing {0} and associates".format(linkname)
        linkbase, ext = os.path.splitext(linkname)
        for exta in self.groupdict.get(ext, [ext]):
            # all related extensions
            try:
                os.remove(linkbase + exta)
            except:
                # if file does not exist ignore the error
                pass

        try:
            os.remove(linkbase + ".o")
        except OSError:
            pass

        try:
            os.remove(linkbase.lower() + ".mod")
        except OSError:
            try:
                os.remove(linkbase.upper() + ".mod")
            except OSError:
                pass

    def is_setup_file(self, fname):
        """return if fname is in NO_CLOBBER_EXCEPTION_LIST"""
        fname = os.path.basename(fname)
        for prefix in setup_globals.NO_CLOBBER_EXCEPTION_LIST:
            if fname[:len(prefix)] == prefix:
                return 1
        return 0

    def clean_setup_dir(self):
        if not os.path.isdir(self.setup_dir):
            os.makedirs(self.setup_dir)

        if gvars.noClobber:
            # dont delete all files
            gvars.out.put('removing files except source and object files (if any)',
                          setup_globals.INFO)
        else:
            gvars.out.put('removing old links in build directory {0}'.format(self.setup_dir),
                          setup_globals.INFO)

        for fname in os.listdir(self.setup_dir):
            # ignore directories
            if os.path.isdir(os.path.join(self.setup_dir, fname)):
                continue

            # what extension?
            if (gvars.noClobber and (os.path.splitext(fname)[1] in self.no_clobber_ext)):
                continue
            else:
                # should we keep this file
                if not self.is_setup_file(fname) or (os.path.splitext(fname)[1] == '.o'):
                    os.remove(os.path.join(self.setup_dir, fname))

    def check_links(self):
        # check to see if we are trying to link in a ABC.F90 as well as an ABC.C, if so
        # inform user of potential problem and quit
        keys = self.links.keys()
        for file in keys:
            base, ext = os.path.splitext(file)
            be = base + ext
            for ext2 in self.groupdict[ext]:
                be2 = base + ext2
                if ext2 != ext and be2 in keys:
                    msg = ('ERROR Checking Links: Both {0} and {1} should not be in '
                           'object directory')
                    msg = msg.format(be, be2)
                    raise SetupError(msg)

    def samefile(self, filea, fileb):
        """The following method should check if the two files are the same by computing
        the adler32 checksum. (It is theoretically possibile but UNLIKELY that two different
        text files produce the same checksum.) However, this version of link_files.py is
        adapted for machines where the zlib module of Python is unavailable, so if we
        cannot use it (as indicated by a NameError exception), let the files appear different."""
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
        if not cksuma:
            return None

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

    def changed_links(self):
        """Find the list of files in the current dir with significant extension
        compare with self.links
        links no longer present and are not derived files get deleted
        for links where target has changed the .o file gets deleted
        target same and has been updated is no problem (make takes care of it)"""
        for fname in os.listdir(gvars.project_setup_dir):
            fname = os.path.join(gvars.project_setup_dir, fname)
            # Dont touch files which should be retained between runs
            # or handled elsewhere in the code
            if self.is_setup_file(fname):
                continue

            # Do not touch derived files  (.o,.mod)
            if os.path.splitext(fname)[1] in self.der_ext:
                continue
    
            # the place from where the file should be copied (for the new run)
            # this is None if file not required for new run
            newtarget = self.links.get(fname, None)
            if newtarget:
                if os.path.islink(fname):
                    # fname is a symlink, so check if it points to new target
                    same = (os.path.abspath(os.readlink(fname)) == newtarget)
                else:
                    same = self.samefile(fname, newtarget)
            else:
                same = None

            # target changed or the link is no longer required
            if not same:
                msg = "current file [{0}] is different from [{0}]. Removing it"
                msg = msg.format(fname, newtarget)
                #print msg
                self.rm_link(fname)
