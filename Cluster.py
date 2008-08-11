#!/usr/bin/python
import sys, os
import tables, numpy

filename = sys.argv[1]

h1 = tables.openFile(filename)
particles = h1.root.Particles.All

p = [particles.read()["x"], particles.read()["y"], particles.read()["z"]]
natoms = h1.root.Properties.natoms
print natoms
List = numpy.array((natoms,1))
