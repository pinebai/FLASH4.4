#!/usr/bin/env python

# This script updates the FLASH3 code repository where it is located

# this script requires python 2.4 or later for the subprocess module
import sys, os
from subprocess import Popen, PIPE

SVN='/opt/CollabNet_Subversion/bin/svn'

# this script assumes it is running inside FLASH3/tools/scripts/
pathToFlash = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))
#print "pathToFlash is ", pathToFlash

# update this copy of Flash
args = [SVN, "update", "-q", pathToFlash]
#args = ["svn", "update",  pathToFlash] #remove the -q argument for testing
p = Popen(args, stderr=PIPE)
exitStatus = p.wait()

if exitStatus != 0:
  print "unable to update"
  print p.stderr.read()
  sys.exit(1)

# else
scriptsToExecute = ["codeCheck.py",
                    "rpDoc.py",
                    "fixRobodocHeaders.py"]
errors = []

for scriptToExecute in scriptsToExecute:
  p = Popen(os.path.join(pathToFlash, "tools/scripts", scriptToExecute), stderr=PIPE)
  exitStatus = p.wait()
  if exitStatus != 0:
    errors.append(p.stderr.read())

if errors:
  # revert the code
  args = [SVN, "revert", "-q", "-R", os.path.join(pathToFlash, "source")]
  p = Popen(args, stderr=PIPE)
  exitStatus = p.wait()
  # print the error and exit
  print "\n".join(errors)
  sys.exit(1)

# else all scripts ran sucessfully - recommit
args = [SVN, "commit", "-m", "automated nightly commit by update.py", pathToFlash]
p = Popen(args, stderr=PIPE)
exitStatus = p.wait()

if exitStatus != 0:
  print "unable to commit"
  print p.stderr.read()
  args = [SVN, "revert", "-q", "-R", os.path.join(pathToFlash, "source")]
  p = Popen(args, stderr=PIPE)
  exitStatus = p.wait()
  sys.exit(1)
