#!/usr/bin/python

import sys, string, os, os.path, popen2

# This is a script that searches a directory structure to find
# and replace strings.  We are specifically using it to change the
# names of our API interface routine names for the Grid Unit, but
# it is generic and can be used anywhere

#This script can be used to change names, do an svn move.  You have to
#choose the appropriate routine in "main" at the bottom of this file.

#usage
# > python replacescript.py Grid_names.txt

#where Grid_names.txt is a file in the format
#oldname <\t> newname  -- with as many lines as you want.

#Grid_getRange    Grid_getBlkIndexLimits
#Grid_stupidName     Grid_betterName

class replace:

    def open_file(self):
        filename = sys.argv[1]
        print filename
        input = open(filename, 'r')
        return input


    def read_input(self, input):
        namelist = []

        
        for line in input.readlines():
            entry = []
            entry = string.split(line)
            namelist.append(entry)

        print namelist    
        return namelist    


    #pass in a name list and logical value, to exeNow =1 or test with a grep
    #command, exeNow=0
    def replace_name(self, namelist, fileExt, exeNow):

        for i in range(len(namelist)):
            oldname = namelist[i][0]
            newname = namelist[i][1]

            #this is the regular expression I use, dan uses the \b to get the beginning of a word
            beginStr = 'find . -name ' + fileExt + ' | xargs perl -pi -e "s/([ ,:\\t!\/])'
            endStr = oldname + '/\\1' + newname + '/"'
            exeStr = beginStr + endStr
            print exeStr
            if exeNow:
                print 'executing now!'
                os.system(exeStr)
            else:
                #change file type from "*.F90" to "Makefile" or whatever you want
                grepStr = 'find . -name ' + fileExt + ' | xargs grep ' + oldname
                os.system(grepStr)


    #pass in a name list and logical value, to exeNow =1 or test with a grep
    #command, exeNow=0 -- This version gets names right at edge of file
    def replace_nameMakefile(self, namelist, exeNow):
        
        for i in range(len(namelist)):
            oldname = namelist[i][0]
            newname = namelist[i][1]
            beginStr = 'find . -name "Makefile" | xargs perl -pi -e "s/^'
            endStr = oldname + '/' + newname + '/"'
            exeStr = beginStr + endStr
            print exeStr
            if exeNow:
                print 'executing now!'
                os.system(exeStr)
            else:
                grepStr = 'find . -name "Makefile" | xargs grep ' + oldname
                os.system(grepStr)

                          


    #pass in a name list and logical value, to exeNow =1 or test with a grep
    #command, exeNow=0.  This routine runs the actual svn move command
    def svn_mv(self, namelist, exeNow):

        for i in range(len(namelist)):
            oldname = namelist[i][0] + ".F90"
            newname = namelist[i][1] + ".F90"

            grepStr = 'find . -name ' + oldname 

            #Noel's popen command that captures std out from child process
            p = popen2.Popen3(grepStr)
            childStdOut = p.fromchild.read()
            p.wait()
            

            pathList = [ x.strip('./') for x in childStdOut.split('\n') ]
            #print pathList

            current_path = os.getcwd()
            for j in range(len(pathList)-1):
                oldfile = os.path.join(current_path, pathList[j])
                head, tail = os.path.split(oldfile)
                newfile = os.path.join(head, newname)
                            
                exeStr= "svn mv --force " + oldfile + " " + newfile
                print exeStr
                
                if exeNow:
                    print 'executing now!'
                    os.system(exeStr)
                else:
                    os.system(grepStr)
               


###########################################################################
#This is the main routine

#Get instance of replace class
r = replace()

#open the replace names text file
fileinput = r.open_file()

#read the file into a python list
namelist = r.read_input(fileinput)

#replace names in all *.F90 files, you can change this
#the last arg can be '1' or '0', 1 says execute now!!!
#0 is for testing to verify that your regular expression is correct.  It will
#execute a grep instead.
r.replace_name(namelist, '"*.F90"', 1)

#basically same as above but slightly different regular expression for Makefiles
#to get beginning of line
#r.replace_nameMakefile(namelist, 1)

#routine will do an svn move to changes names of files in repos
#r.svn_mv(namelist, 1)




# search directory structure and replace string
#old stuff


# noel's original line
# find . -name "*.F90" | xargs perl -pi -e "s/([ ,:(])ihi/\1io_ihi/"


# noel's original line
# find . -name "*.F90" | xargs perl -pi -e "s/^ihi/io_ihi/"




