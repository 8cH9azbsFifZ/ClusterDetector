#!/usr/bin/python
import sys, os
import tables, numpy
import cluster
import numpy

import hdf

hh=hdf.hdf("dump.test.h5")

c=cluster

pp=hh.particles.All.readWhere("y>30")
xx = pp["x"]
yy = pp["y"]
zz = pp["z"]
x=numpy.array(xx)
y=numpy.array(yy)
z=numpy.array(zz)

cc=cluster.clusterdetector(x,y,z,3.615*2.)
print cc
