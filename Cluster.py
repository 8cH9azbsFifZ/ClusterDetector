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
      print self.natoms,"particles"

   def DetermineExtrema(self):
      print "Determining extrema"
      self.xlo = self.x.min()
      self.xhi = self.x.max()
      self.ylo = self.y.min()
      self.yhi = self.y.max()
      self.zlo = self.z.min()
      self.zhi = self.z.max()
      self.lx = self.xhi-self.xlo
      self.ly = self.yhi-self.ylo
      self.lz = self.zhi-self.zlo

   def InitList(self):
      print "Init list"
      self.List = numpy.array([i for i in range(0, self.natoms) ])

   def ClustereStraight(self):
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

   def CoordinateToCell(self, x, y, z):
      return int( numpy.floor (x*self.CellLengthI[X]) + numpy.floor (y*self.CellLengthI[Y])*self.NCellSide[X] + numpy.floor (z*self.CellLengthI[Z])*self.NCellSide[X]*self.NCellSide[Y]   )


   def BuildLinkedlist(self):
      CellSize = self.RcNeighbor
      
      ncellx = int (numpy.ceil (self.lx / CellSize))
      ncelly = int (numpy.ceil (self.ly / CellSize))
      ncellz = int (numpy.ceil (self.lz / CellSize))
      ncell = ncellx+ncelly+ncellz

      CellLengthx = self.lx / ncellx
      CellLengthy = self.ly / ncelly
      CellLengthz = self.lz / ncellz

      CellLengthIx = 1./CellLengthx
      CellLengthIy = 1./CellLengthy
      CellLengthIz = 1./CellLengthz

      for i in range (0, self.natoms):

         LinkedListInsert (i, CellList + CoordinateToCell (x[i], y[i], z[i]) );


   def ClustereVerletlist(self):
      print "Clustere"
      for i in range (0, self.natoms):
         if i != self.List[i]:
            continue

         j = i

         while True:
            # sum overall neighbor cells
            for mx in [-1, 0, 1]:
               x = self.x[j] + mx*RcCluster
               if x < self.xlo or x > self.xhi:
                  continue
               for my in [-1, 0, 1]:
                  y = self.y[j] + my*RcCluster
                  if y < self.ylo or y > self.yhi:
                     continue
                  for mz in [-1, 0, 1]:
                     z = self.z[j] + mz*RcCluster
                     if z < self.zlo or z > self.zhi:
                        continue

                     # sum overall particles in cell

            j = self.List[j]

            if i == j:
               break

   def __init__(self, RcCluster, filename):
      self.RcCluster = RcCluster
      self.RcClusterSq = RcCluster*RcCluster

      self.ReadFile (filename)
      self.DetermineExtrema ()
      self.InitList ()


c = Cluster(3.615, filename)
