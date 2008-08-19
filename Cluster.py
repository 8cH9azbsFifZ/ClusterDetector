#!/usr/bin/python
import sys, os
import tables, numpy
import cluster
import numpy

import hdf

hh=hdf.hdf("dump.test.h5")

c=cluster

x = hh.particles.All.read()["x"]
y = hh.particles.All.read()["y"]
z = hh.particles.All.read()["z"]

xx=numpy.array(x)
yy=numpy.array(y)
zz=numpy.array(z)

cc=cluster.clusterdetector(xx,yy,zz,3.615*2.)
print cc
