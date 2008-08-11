#!/usr/bin/python
import sys, os
import tables, numpy

filename = sys.argv[1]



def SqDist(a,b):
   return (a-b)**2

class Cluster:
   def ReadFile(self, filename):
      print "Read file", filename
      self.filename = filename
      self.h1 = tables.openFile (self.filename)
      particles = self.h1.root.Particles.All
      [self.x, self.y, self.z] = [particles.read()["x"], particles.read()["y"], particles.read()["z"]]
      self.natoms = numpy.size (self.x)

   def InitList(self):
      print "Init list"
      self.List = numpy.array ((self.natoms,1))
      for i in range(0,natoms):
         self.List[i] = i

   def Clustere():
      print "Clustere"
      for i in range(0, self.natoms-1):
         if i != self.List[i]:
            continue

         j = i
         Xj = self.x[j]
         Yj = self.y[j]
         Zj = self.z[j]

         while True:
            for k in  range(i+1, self.natoms):
               Lk = self.List[k]
               if k == Lk:
                  RSq = SqDist(Xj, self.x[k]) + SqDist(Yj, self.y[k]) + SqDist(Zj, self.z[k])
                  if RSq <= self.RcClusterSq:
                     self.List[k] = self.List[j]
                     self.List[j] = Lk
            
            j = self.List[j]
            Xj = self.x[j]
            Yj = self.y[j]
            Zj = self.z[j]
            
            if i == j:
               break


c = Cluster()
c.ReadFile (filename)
c.InitList()

