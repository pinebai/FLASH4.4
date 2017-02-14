#!/usr/bin/env python
import re, shutil, glob, os, sys

class UnitMakefile:
    def __init__(self):
        self.objList = []
        self.dependencies = {}
        #Preprocessor sections are only copied and pasted.  This should
        #be OK for now as we currently only have MODUPPERCASE sections that
        #deal with module file dependencies.  We may need increased
        #sophistication in future.
        self.preProcess = ""
        self.unit = ""

    #Note that there are many limitations to readFromFile function.  It does
    #not handle all of the permutations of make language, especially line
    #continuation characters within my pattern strings.  Fortunately if
    #it messes up it will only create a broken Makefile and so the error
    #is easy to find.
    def readFromFile(self,unit,name):
        self.unit = unit

        #http://code.google.com/edu/languages/google-python-class/regular-expressions.html
        #pat2 is the same as pat4 except pat4 will detect file extensions
        #beyond .o which is necessary for dependencies.

        #Used for build list.
        #pat1 allows us to detect "Unit +=" or "Unit =" build lists.
        pat1 = re.compile("^\s*%s\s*\+?=.*" % unit, re.IGNORECASE)
        pat2 = re.compile("\w+\.o", re.IGNORECASE)

        #Used for dependencies.
        #Note that pat4 allows us to depend on e.g. a .o or a .h or a .fh file.
        pat3 = re.compile("^\s*\w+\.o\s*:", re.IGNORECASE)
        pat4 = re.compile("\w+\.\w+", re.IGNORECASE)

        #Used for preprocessor sections.
        #pat5 has extra parenthesis and so we get a list of tuples.
        pat5 = re.compile("(^\s*ifdef\s+(\w+))", re.IGNORECASE)
        pat6 = re.compile("^\s*endif\s*$", re.IGNORECASE)


        if not (os.path.isfile(name)):
            print "File ", name, " does not exist"
            sys.exit(1)

        file = open(name, "r")
        lines = file.read().split("\n")
        file.close()

        i = 0
        make_lines = []
        depend_lines = []

        while i < len(lines):

            m1 = pat1.findall(lines[i])
            m3 = pat3.findall(lines[i])
            m5 = pat5.findall(lines[i])
            line = ""

            if (m1):
                if lines[i].endswith("\\"):
                    while lines[i].endswith("\\"):
                        line += lines[i][:-1]
                        i+=1
                    line += lines[i]
                else:
                    line = lines[i]
                make_lines.append(line)


            elif (m3):
                if lines[i].endswith("\\"):
                    while lines[i].endswith("\\"):
                        line += lines[i][:-1]
                        i+=1
                    line += lines[i]
                else:
                    line = lines[i]
                depend_lines.append(line)

            elif (m5):
                if (m5[0][1] <> "MODUPPERCASE"):
                    print "WARNING - Macro is not MODUPPERCASE"
                while i < len(lines):                        
                    self.preProcess += lines[i] + '\n'
                    m = pat6.findall(lines[i])
                    if (m):
                        break
                    i+=1

            i+=1

        self.objList = [m2 for line in make_lines for m2 in pat2.findall(line)]
        #This is here because a Makefile may contain multiple references to the
        #same object file.
        self.objList = list(set(self.objList))

        self.dependencies = {}
        for line in depend_lines:
            m4 = pat4.findall(line)
            key = m4[0] #First file name in list of all file names
            values = m4[1:] #Remaining file name list

            if self.dependencies.has_key(key):
                prior_values = self.dependencies[key]
                values.extend(prior_values)
            self.dependencies[key] = values

    def writeToFile(self,name):
        makefile = open(name, 'w')
        makefile.write("# Object list\n" + self.unit + " +=")
        for obj in self.objList:
            makefile.write(" " + obj)
        makefile.write("\n")

        makefile.write("\n# Dependency list\n")
        for lhs in self.dependencies.keys():
            rhs = self.dependencies[lhs]
            makefile.write(lhs + " :")
            for x in rhs:
                makefile.write(" " + x)
            makefile.write("\n")

        makefile.write("\n# Preprocessor section\n")
        makefile.write(self.preProcess)
        makefile.close()

    def removeObjects(self,remList):
        #Amend the Makefile object list.
        origList = list(self.objList)
        self.objList = [i for i in origList if i not in remList]

        #Remove dependencies where we have removed the LHS object.
        keys = self.dependencies.keys()
        for obj in remList:
            if obj in keys:
                print "Delete target", obj, "which has dependencies", \
                    self.dependencies[obj]
                del self.dependencies[obj]

        #Amend dependencies where we have removed the RHS object.
        keys = self.dependencies.keys()
        for obj in keys:
            deps = self.dependencies[obj]
            remDeps = [x for x in deps if x in remList]
            if remDeps:
                print "Remove dependencies", remDeps, \
                    "for target", obj, "which had dependencies", \
                    self.dependencies[obj]
                newDeps = list(set(deps) - set(remDeps))
                if newDeps:
                    self.dependencies[obj] = newDeps
                else:
                    del self.dependencies[obj]

    def getObjList(self):
        return self.objList

    def getDependencies(self):
        return self.dependencies

    def printInfo(self):
        print "Object List:\n", self.objList
        print ""
        print "Dependency List:\n", self.dependencies
        print ""
        print "Preprocessor Lines:\n", self.preProcess


if __name__=="__main__":
    makefileName = "Makefile.Grid"
    unitName = "Grid"

    unitMakefile = UnitMakefile()
    unitMakefile.readFromFile(unitName,makefileName)
    unitMakefile.printInfo()
    unitMakefile.writeToFile(makefileName + ".NEW")
