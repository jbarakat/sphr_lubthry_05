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
outdir   = workdir + '/../output/go/Ca1em01/Bo1ep00/Ma1ep01/dr5em02/dt1em05'
#outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em05'

outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma1ep01/tstop1-2/dr5em02/dt1em05'
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
time    = ax1.text(0.85, 0.05, '', transform=ax1.transAxes)

for ax in (ax1,ax2, ax3, ax4):
	ax.set_xlim( 0, 5)
	ax.set_xlabel('$\overline{\sigma}$')
ax1.set_ylabel('$\overline{h}$')
ax1.set_xlim( 0, 1)
ax1.set_ylim( 0, 0.3)
ax2.set_ylabel('$\overline{{\Gamma}}$')
ax2.set_ylim(-0.5, 1)
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
	time  .set_text('')
	return line1a, line1b, line2, line3, line4, time

# animation function, to be called sequentially
def animate(i) :
	t   =  tdata[0,i]
	x   =  rdata[:,i]
	y1a =  hdata[:,i]
	y1b = -fdata[:,i]
	y2  =  gdata[:,i] + 1
	y3  =  qdata[:,i]
	y4  =  vdata[:,i]

	line1a.set_data(x, y1a)
	line1b.set_data(x, y1b)
	line2 .set_data(x, y2 )
	line3 .set_data(x, y3 )
	line4 .set_data(x, y4 )
	time  .set_text("t = " + str(t))
	return line1a, line1b, line2, line3, line4, time

# call the animator. blit=True means only re-draw the parts that have changed
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=N1, interval=30, repeat=True, blit=True)

## format and save movie
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
#anim.save('im.mp4', writer=writer)
                               #fargs=(rdata, hdata, -fdata, 1+gdata, qdata, vdata,
															 #       line1a, line1b, line2, line3, line4, time),

plt.show()
