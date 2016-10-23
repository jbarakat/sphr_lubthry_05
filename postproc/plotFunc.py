import numpy as np
from numpy import column_stack
from matplotlib import pyplot as plt 
from matplotlib.pyplot import *
from math import pi
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
outdir   = workdir + '/../output'
tfile    = '/t.txt'
rfile    = '/r.txt'
hfile    = '/h.txt'
gfile    = '/g.txt'
ffile    = '/f.txt'
pfile    = '/p.txt'
qfile    = '/q.txt'
vfile    = '/vs.txt'

# load data
tdata    = np.loadtxt(outdir + tfile, unpack=True,skiprows=0)
rdata    = np.loadtxt(outdir + rfile, unpack=True,skiprows=0)
hdata    = np.loadtxt(outdir + hfile, unpack=True,skiprows=0)
gdata    = np.loadtxt(outdir + gfile, unpack=True,skiprows=0)
fdata    = np.loadtxt(outdir + ffile, unpack=True,skiprows=0)
pdata    = np.loadtxt(outdir + pfile, unpack=True,skiprows=0)
qdata    = np.loadtxt(outdir + qfile, unpack=True,skiprows=0)
vdata    = np.loadtxt(outdir + vfile, unpack=True,skiprows=0)

# number of space and time points
J1 = len(tdata[:,0])
N1 = len(tdata[0,:])
J  = J1 - 1
N  = N1 - 1

pl = 'x'
ydata = hdata

if pl == 'x' :
	for i in range(N1) :
		if i % 1 == 0 :
			x = rdata[:,i]
			y = ydata[:,i]
			plot(x,y,'-')
		
		#y = -fdata[:,i]
		#plot(x,y,'k--')
	#xlim([0,5])
	#ylim([0,max(y)+0.01])
	show()

if pl == 't' :
	for i in range(J1) :
		if i % 1 == 0 :
		#if i < 20 :
			x = tdata[i,:]
			y = ydata[i,:]
			plot(x,y,'-')
	show()

