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
outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma1ep03/dr5em02/dt1em05/'
#outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em05/'
outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em04/'

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
fig = plt.figure(figsize=(20,10))
#ax = plt.axes(xlim=(0, 5), ylim=(-2,2))
#line, = ax.plot([], [], lw=2)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
line1a, = ax1.plot([], [], 'b-')
line1b, = ax1.plot([], [], 'k-')
line2,  = ax2.plot([], [], 'g-')
line3,  = ax3.plot([], [], 'b-')
line4,  = ax4.plot([], [], 'g-')

for ax in (ax1,ax2, ax3, ax4):
	ax.set_xlim( 0, 20)
	ax.set_xlabel('$\overline{\sigma}$')
ax1.set_ylabel('$\overline{h}$')
ax1.set_ylim(-1.5, 1.5)
ax2.set_ylabel('$\overline{{\Gamma}}$')
ax2.set_ylim( 0, 1)
ax3.set_ylabel('$\overline{q}$')
ax3.set_ylim( -0.5, 10)
ax4.set_ylabel('$\overline{v}_s$')
ax4.set_ylim( -0.5, 1)


# initialization function: plot the background of each frame
def init() :
	line1a.set_data([], [])
	line1b.set_data([], [])
	line2 .set_data([], [])
	line3 .set_data([], [])
	line4 .set_data([], [])
	return line1a, line1b, line2, line3, line4

# animation function, to be called sequentially
def animate(i,xdata,y1adata,y1bdata,y2data,y3data,y4data,
              line1a,line1b,line2,line3,line4) :
	print "ts = " + str(i)
	x   = xdata [:,i]
	y1a = y1adata[:,i]
	y1b = y1bdata[:,i]
	y2  = y2data [:,i]
	y3  = y3data [:,i]
	y4  = y4data [:,i]
	line1a.set_data(x, y1a)
	line1b.set_data(x, y1b)
	line2 .set_data(x, y2 )
	line3 .set_data(x, y3 )
	line4 .set_data(x, y4 )
	return line1a, line1b, line2, line3, line4

# call the animator. blit=True means only re-draw the parts that have changed
anim = animation.FuncAnimation(fig, animate,
                               fargs=(rdata, hdata, -fdata, 1+gdata, qdata, vdata,
															        line1a, line1b, line2, line3, line4),
                               init_func=init, frames=N1, interval=30, repeat=True, blit=True)

## format and save movie
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
#anim.save('im.mp4', writer=writer)

plt.show()
