#!/usr/bin/python -i

# execute this file with two arguments pointing to drift tuple logs.  the two
# logs will be loaded into relations 'a' and 'b' respectively and leave you with
# an interactive python shell.

import sys

# requires Dee: http://www.quicksort.co.uk/DeeDoc.html
from Dee import *

# load a file of tuples and return as a relation
def read_tups(path):
  with open(path) as f:
    def impt(x):
      t = eval('(' + x + ')')
      return (t[0],t[1],t[2],t[3],bool(t[4]),t[5],long(t[6]))
    tups = [impt(x) for x in f]
  return Relation(['inst','step','src','blk','leaf','unk','val'], tups)

# a relation where for each (step,blk,unk) it gives the maximum inst
def latest(r):
  r1 = r(['inst','step','blk','unk'])
  return SUMMARIZE(r1, r1(['step','blk','unk']), {'inst':(MAX,lambda t:t.inst)})

# find everything in a that doesnt have a corresponding (blk,unk,val) in b
def diff(a,b):
  return a - (a & b(['blk','unk','val']))

# notation for generating a singleton relation, ex gen(a=1,b=2)
def gen(**a):
  return GENERATE(a)

a = read_tups(sys.argv[1])
b = read_tups(sys.argv[2])

#print diff(latest(a),latest(b))
