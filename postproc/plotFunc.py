import numpy as np
from numpy import column_stack
from matplotlib import pyplot as plt 
from matplotlib.pyplot import *
from math import pi, log
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
outdir   = workdir + '/../output/go/Ca1em03/Bo1ep00/Ma0ep00/dr5em02/dt1em05/'
outdir   = workdir + '/../output/go/Ca1em03/Bo1ep00/Ma1ep00/dr5em02/dt1em05/'
#outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em02/'

# filenames
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

pl = 'x'
ydata = vdata

# number of space and time points
J1 = len(tdata[:,0])
N1 = len(tdata[0,:])
J  = J1 - 1
N  = N1 - 1


if pl == 'x' :
	for i in range(N1) :
		if i % 5 == 0 :
		#if i < 100 :
			x = rdata[:,i]
			y = ydata[:,i]
			plot(x,y,'-')
		
		#y = -fdata[:,i]
		#plot(x,y,'k--')
	#xlim([0,0.05])
	#ylim([0.2,0.21])
	
	#xlim([0,5])
	#ylim([0,max(y)+0.01])
	show()

if pl == 't' :
	for i in range(J1) :
		#if i % 1 == 0 :
		if i == 0 :
			x = tdata[i,:]
			y = ydata[i,:]
			for n in range(N1) :
				y[n] = hdata[i,n] + fdata[i,n]
			#loglog(x,y,'-')
			plot(x,y,'-')
	xlim([1.2,50])
	ylim([0,0.01])
	show()

