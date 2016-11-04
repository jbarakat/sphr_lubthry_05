import numpy as np
from numpy import column_stack
from matplotlib import pyplot as plt 
#from matplotlib.pyplot import *
from matplotlib import animation
from math import pi, log
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
workdir  = workdir + '/..'
outdir   = workdir + '/../output/go/Ca1em01/Bo1ep00/Ma0ep00/dr5em02/dt1em05'
outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em06'
#outdir   = workdir + '/../output/go/Ca1em03/Bo1ep00/Ma0ep00/dr5em02/dt1em06'
#outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em05'

#outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em04'
#outdir   = workdir + '/../output/stop/Ca1em01/Bo1ep00/Ma1ep02/tstop1-2/dr5em02/dt1em4'

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

# number of space and time points
J1 = len(tdata[:,0])
N1 = len(tdata[0,:])
J  = J1 - 1
N  = N1 - 1

# set up the figure, the axis, and the plot element to be animated
fig = plt.figure(figsize=(15,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax2 = fig.add_subplot(133)
line1a, = ax1.plot([], [], 'b-')
line1b, = ax1.plot([], [], 'k-')
line2a, = ax2.plot([], [], 'b-')
line2b, = ax2.plot([], [], 'k-')
line3a, = ax3.plot([], [], 'b-')
line3b, = ax3.plot([], [], 'k-')
time1   = ax1.text(0.85, 0.05, '', transform=ax1.transAxes)
time2   = ax2.text(0.85, 0.05, '', transform=ax2.transAxes)
time3   = ax3.text(0.85, 0.05, '', transform=ax3.transAxes)

for ax in (ax1,ax2):
	ax.set_xlim( 0, 4)
	ax.set_xlabel('$\overline{\sigma}$')
	ax.set_ylabel('$\overline{h}$')
	ax.set_ylim(-1.5, 1.5)

# initialization function: plot the background of each frame
def init() :
	line1a.set_data([], [])
	line1b.set_data([], [])
	line2 .set_data([], [])
	time1 .set_text('')
	time2 .set_text('')
	return line1a, line1b, line2, time1, time2

# animation function, to be called sequentially
def animate(i) :
	t   =  tdata[0,i]
	x   =  rdata[:,i]
	y1a =  hdata[:,i]
	y1b = -fdata[:,i]
	y2  =  gdata[:,i] + 1

	line1a.set_data(x, y1a)
	line1b.set_data(x, y1b)
	line2 .set_data(x, y2 )

	time1 .set_text("$\overline{t}$ = " + str(t))
	time2 .set_text("$\overline{t}$ = " + str(t))

	print "t = " + str(t)
	return line1a, line1b, line2, time1, time2

# call the animator. blit=True means only re-draw the parts that have changed
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=N1, interval=45, repeat=True, blit=True)

# format and save movie
fileName = 'hg_Ca1em01_Bo1ep00_Ma0ep00_dr5em02_dt1em05.mp4'
fileName = 'hg_Ca1em02_Bo1ep00_Ma0ep00_dr5em02_dt1em06.mp4'
Writer = animation.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
anim.save(fileName, writer=writer)

plt.show()
