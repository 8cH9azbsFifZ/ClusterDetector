#!/usr/bin/python
import sys, os
import tables, numpy
import cluster

c=cluster

x = [1,2,3,1,2,3,4,10]
y = [2,3,4,2,3,4,5,3]
z = [1,1,2,5,2,5,7,11]

print c.clusterdetector(x,y,z,3.615*2.)

