#!/usr/bin/env python

import re, sys

reobj = re.compile(r"[a-z]",re.I)
partline = re.compile(r"\\\s*$",re.I)

def process(filename):
  global reobj
  fd = file(filename)
  remove = 0
  for line in fd.readlines():
      if line.find(":") >= 0: remove = 1
      if remove == 0: print line,
      # Is this line being continued
      ob = partline.search(line)
      if ob is None: remove = 0

if __name__=="__main__":
   if len(sys.argv) < 2:
      print >> sys.stderr, "Usage: cleanup.py <Makefile> > NewMakefile"
   else: 
      process(sys.argv[1])

