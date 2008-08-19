#!/usr/bin/python
import sys, os
import tables, numpy
import cluster

import hdf

hh=hdf.hdf("dump.test.h5")

c=cluster

x = hh.particles.All.read()["x"]
y = hh.particles.All.read()["y"]
z = hh.particles.All.read()["z"]

print c.clusterdetector(x,y,z,3.615*2.)

