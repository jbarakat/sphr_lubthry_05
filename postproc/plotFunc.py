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
outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em05/'

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

pl = 't'
ydata = hdata

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

if pl == 't' :
	for i in range(J1) :
		#if i % 1 == 0 :
		if i == 0 : # centerline
			x = tdata[i,:]
			y = ydata[i,:]
			for n in range(N1) :
				y[n] = hdata[i,n] + fdata[i,n]
				#y[n] = gdata[i,n]+1
			#loglog(x,y,'-')
			loglog(x,y,'-')

			# calculate slope of log log plot
			logx = []
			logy = []
			for n in range(N1) :
				if y[n] > 0 and x[n] > 0:
					logx.append(log(x[n]))
					logy.append(log(y[n]))
				else :
					logx.append(log(x[1]))
					logy.append(log(y[1]))
			#plot(logx,logy,'-')

			slope = (logy[-1] - logy[-2])/(logx[-1] - logx[-2])
			print slope

	#xlim([1.2,50])
	#ylim([0,0.01])
	#ylim([0,0.001])
	xlabel('$\overline{t}$')
	ylabel('$\overline{h}_1 + \overline{h}_2$')
	#ylabel('$\overline{\Gamma}$')
show()

