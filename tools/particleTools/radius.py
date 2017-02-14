 #!usr/bin/python

from math import *
import thread, time, string, sys





#This program takes 2 particle dump files and finds the distance between
#particles.  It is really only applicable to the Orbit problem

#usage > python radius.py test_ParticlesDump1_0000 test_ParticlesDump2_0000

#output is to stdout with data in 4 colums in the format
# <nstep> <dist between 2 particles in 1st file> <dist between 2 particles in 2nd file> <diff>

class comp:

  #main driver for this program
  def run(self):
    print "inside run"
    input1, input2 = comp.open_files(self)

    data1 = comp.read_file(self, input1)
    data2 = comp.read_file(self, input2)

    dist1 = comp.analyze(self, data1)
    dist2 = comp.analyze(self, data2)

    comp.compare(self, dist1, dist2)
    
    comp.close_file(self, input1)
    comp.close_file(self, input2)


  #open the files  
  def open_files(self):
    filename1 = sys.argv[1]
    print filename1
    input1 = open(filename1, 'r')

    filename2 = sys.argv[2]
    print filename2
    input2 = open(filename2, 'r')

    return input1, input2




  def close_file(self, input):
    input.close()

    
  #read data file and return a list 
  def read_file(self, input):
    part_data = []

    for line in input.readlines():
      part = line.split()
      part_data.append(part)

    return part_data


  #compute the distance between the 2 particles
  def analyze(self, data):

    dist = []

    for i in range(0,len(data),2):
      x1 = float(data[i][2])
      y1 = float(data[i][3])
      x2 = float(data[i+1][2])
      y2 = float(data[i+1][3])

      distance = sqrt((x2-x1)**2 + (y2-y1)**2)
      dist.append(distance)
      i = i+2
      
    return dist


  #compare 2 different files
  def compare(self, dist1, dist2):
    
    for i in range(len(dist1) -1):
      print i, dist1[i], dist2[i], dist1[i] - dist2[i]



########################################################################

#print arguments
print sys.argv

#instantiate class    
c = comp()

#run program
c.run()
