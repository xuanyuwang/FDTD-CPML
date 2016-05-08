# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

dz = 0.015
space = 16 
time =  150
data = np.loadtxt('ex.txt')

fig = plt.figure()
ax = plt.axes(xlim = (0,1), ylim = (-0.1,1))
line, = ax.plot([],[],lw = 2)

def init():
	line.set_data([],[])
	return line,

def animate(i):
    x = np.arange(0 * dz, (space + 1) * dz, dz)
    y = data[i]
    line.set_data(x,y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=time, interval=50, blit=True)

#anim.save('ba.mp4', fps = 30, extra_args=['-vcodec', 'libx264'])

plt.show()