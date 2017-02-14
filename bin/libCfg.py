
## Class to parse CONFIG files for library

__all__ = [ "FlashLib"]

##############################

import string, re, UserDict, types, os

import globals
from globals import * # GVars and SetupError
from utils import *

class FlashLib(UserDict.UserDict):
    """Encapsulates library information as expressed in Config files.
    Data is accessed through dictionary methods, ex:
  
    FlashLib('pfft')['LIBRARY']

    Currently the only supported keyword in Config files in lib directory is LIBRARY keyword

    NOTE: we need to be in "lib" directory when instantiating this
          All the code are copied from FlashUnit (this code is a subset of that of FlashUnit)
    """

    def __init__(self, pathName, ignoreConfig=0):
        UserDict.UserDict.__init__(self)

        self.COMMENT = '#'
        self.QUOTE   = '"'
        self.FILEBASE = 'Config'

        self.regexps = {}
        self.parsers = []
        pathname = pathName.lower()
        
        # This is just a clever way to list the methods in the
        # class. The dir(self.__class__) lists all the attributes
        # like __cmp__ , __init___, __setitem__, __len__ etc and
        # then also the user defined methods, in this case,
        # getParent, initParser, match, parseDEFAULT, parseEXCLUSIVE etc.
        for name in dir(self.__class__):
            if not re.compile(r'parse[A-Z]+$').match(name): continue
            if type(getattr(self, name))==types.MethodType:
                self.parsers.append(name)
                try: getattr(self, name)('') #initialize regexps, dictionary
                except SetupError: pass
        
        self.update( {'LIBRARY':{}, 'TYPE':'EXTERNAL'})

        if not os.path.isdir(pathname): # given name does not make sense
           GVars.out.put('Directory %s not found. Assuming external library' % os.path.join(os.getcwd(),pathname), globals.DEBUG)
           ignoreConfig = 1
       
        self.name = os.path.normpath(pathname) #something like 'pfft'

        if (not ignoreConfig) and \
          os.path.isfile(os.path.join(self.name, self.FILEBASE)):
            self.filename = os.path.join(self.name, self.FILEBASE)
            self.parse()
        else:
            self.filename = ''

    def __cmp__(self, other):
        """Alphabetical comparison on unit names (like 'source/io/amr')"""
        if type(other)==types.StringType:
            return cmp(self.name, other)
        else:
            return cmp(self.name, other.name)

    def __repr__(self):
        return self.name
   
    def match(self, keyword, line):
        match = self.regexps[keyword].match(line)
        if not match:
            raise SetupError('input doesn\'t match regular expression "%s"'%\
                             self.regexps[keyword].pattern)
        return match

    def parse(self):
        lineno = 0

        for line in open(self.filename).readlines():
            lineno += 1
            rawline = line
            if rawline and rawline[-1]=='\n': rawline = rawline[:-1]

            line=stripComments(line, self.COMMENT,self.QUOTE)
            line=string.strip(line)
            if not line: continue

            pkeyword = "parse%s" % string.split(line)[0]
            if pkeyword not in self.parsers:
                raise SetupError('Unknown keyword: file %s, line '\
                                 '%d\n%s'%(self.filename, lineno, rawline))
            try:
                getattr(self,pkeyword)(line)
            except SetupError, msg:
                raise SetupError('Bad syntax: file %s, line %d:\n%s\n\n%s' % \
                                 (os.path.join("lib",self.filename), lineno, rawline, str(msg)))

    def initParser(self, keyword, initvalue, regexp=None):
        if self.has_key(keyword): return
        self[keyword]=initvalue
        if regexp:
            self.regexps[keyword]=re.compile(regexp)

    def parseLIBRARY(self, line):
        self.initParser('LIBRARY', {}, 'LIBRARY\s+(\S+)\s*(.*)$')
        libmatch = self.match('LIBRARY',line)
        libname = libmatch.group(1).lower()
        libargs = string.join(libmatch.group(2).split()) # trims and removes multiple spaces
        self['LIBRARY'][libname] = libargs

    def parseTYPE(self,line):
        self.initParser('TYPE', "EXTERNAL",'TYPE\s+(INTERNAL|EXTERNAL)\s*$')
        self['TYPE']= self.match('TYPE',line).group(1)

