#!/usr/bin/python
import sys, os
import tables, numpy

import cluster




#filename = sys.argv[1]



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

   def ShiftPositions(self):
      print "Shifting positions"
      if self.xlo < 0:
         self.xhi -= self.xlo
         for i in range(0, self.natoms):
            self.x[i] -= self.xlo
         self.xlo = 0.

      if self.ylo < 0:
         self.yhi -= self.ylo
         for i in range(0, self.natoms):
            self.y[i] -= self.ylo
         self.ylo = 0.

      if self.zlo < 0:
         self.zhi -= self.zlo
         for i in range(0, self.natoms):
            self.z[i] -= self.zlo
         self.zlo = 0.

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
      return int( min(numpy.floor (x*self.CellLengthIx),self.ncellx) + min(numpy.floor (y*self.CellLengthIy)*self.ncellx, self.ncelly) + min(numpy.floor (z*self.CellLengthIz)*self.ncellx*self.ncelly,self.ncellz) )


   def BuildLinkedlist(self):
      print "Build linked list"
      CellSize = self.RcCluster
      
      self.ncellx = int (numpy.ceil (self.lx / CellSize))
      self.ncelly = int (numpy.ceil (self.ly / CellSize))
      self.ncellz = int (numpy.ceil (self.lz / CellSize))

      self.ncell = self.ncellx * self.ncelly * self.ncellz
      print "ncells:",self.ncellx,self.ncelly,self.ncellz,self.ncell

      list = []
      self.LinkedList = [list for i in range(0, self.ncell+1) ]

      self.CellLengthx = self.lx / self.ncellx
      self.CellLengthy = self.ly / self.ncelly
      self.CellLengthz = self.lz / self.ncellz

      self.CellLengthIx = 1./self.CellLengthx
      self.CellLengthIy = 1./self.CellLengthy
      self.CellLengthIz = 1./self.CellLengthz
      print "cell length:",self.CellLengthx,self.CellLengthy,self.CellLengthz,self.CellLengthIx,self.CellLengthIy,self.CellLengthIz

      for i in range (0, self.natoms):
         if i % 10000 == 0:
            print i,
#         print i,self.x[i],self.y[i],self.z[i],self.CoordinateToCell (self.x[i], self.y[i], self.z[i])
         self.LinkedList[self.CoordinateToCell (self.x[i], self.y[i], self.z[i])].append (i) # = [i, self.LinkedList[self.CoordinateToCell (self.x[i], self.y[i], self.z[i])]]


   def ClustereVerletlist(self):
      print "Clustere"
      for i in range (0, self.natoms):
         if i != self.List[i]:
            continue

         j = i

         while True:
            # sum overall neighbor cells
            for mx in [-1, 0, 1]:
               x = self.x[j] + mx*self.RcCluster
               if x < self.xlo or x > self.xhi:
                  continue
               for my in [-1, 0, 1]:
                  y = self.y[j] + my*self.RcCluster
                  if y < self.ylo or y > self.yhi:
                     continue
                  for mz in [-1, 0, 1]:
                     z = self.z[j] + mz*self.RcCluster
                     if z < self.zlo or z > self.zhi:
                        continue

                     # sum overall particles in cell
                     for k in self.LinkedList[self.CoordinateToCell (self.x[j], self.y[j], self.z[j])]:
                        if k % 10000 == 0:
                           print k,
                           sys.stdout.flush()
                        RSq = SqDist(self.x[j], self.x[k]) + SqDist(self.y[j], self.y[k]) + SqDist(self.z[j], self.z[k])
                        if RSq <= self.RcClusterSq:
                            Lk = self.List[k]
                            self.List[k] = self.List[j]
                            self.List[j] = Lk

            j = self.List[j]

            if i == j:
               break

   def Spectrum(self):
      self.Spectrum = numpy.array([0 for i in range(0, self.natoms+1) ])
      for i in range (0, self.natoms):
         if self.List[i] > 0:
            size = 1
            j = self.List[i]
            self.List[i] = -1
            while True:
               size += 1
               k = self.List[j]
               self.List[j] = -1
               j = k
               if j == i:
                  break
            self.Spectrum[size] += 1
      for i in range (0, self.natoms+1):
         if self.Spectrum > 0:
            print i, self.Spectrum


   def __init__(self, RcCluster, filename):
      self.RcCluster = RcCluster
      self.RcClusterSq = RcCluster*RcCluster

      self.ReadFile (filename)
      self.DetermineExtrema ()
      self.ShiftPositions ()
      self.InitList ()
      self.BuildLinkedlist ()
      self.ClustereVerletlist()
      self.Spectrum ()


#c = Cluster(3.615, filename)
