#!/usr/bin/python
import sys, os
import tables, numpy

filename = sys.argv[1]

h1 = tables.openFile(filename)
particles = h1.root.Particles.All

[x, y, z] = [particles.read()["x"], particles.read()["y"], particles.read()["z"]]
natoms = numpy.size (x)
print natoms

def InitList():
   List = numpy.array((natoms,1))
   for i in range(0,natoms):
      List[i] = i

def Clustere():
   for i in range(0,natoms-1):
      if i != List[i]:
         continue
      j = i
      Xj = x[j]
      Yj = y[j]
      Zj = z[j]

