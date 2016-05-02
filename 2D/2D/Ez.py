# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 21:08:46 2016

@author: obser
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import time

'''
import data
'''
dz = 0.015
space =16
ez_space = space + 1
time = 150
data = np.loadtxt('Ez.txt')
groups = np.ndarray(shape=(time, ez_space, ez_space))

for i in range(0, time):
    groups[i] = data[i * ez_space:(i + 1) * ez_space]

'''
set plot data
'''
x = np.arange(0 * dz, ez_space * dz, dz)
y = np.ndarray(shape=(ez_space))
y.fill(0)

X = np.ndarray(shape=(ez_space, ez_space))
Y = np.ndarray(shape=(ez_space, ez_space))
Z = np.ndarray(shape=(ez_space, ez_space))

for i in range(0, ez_space):
    X[i] = x
    Y[i] = y + i * dz

fig = plt.figure()
ax = Axes3D(fig)
ax.set_autoscale_on(True)
ax.set_zlim(-1,1)


'''
plot
'''
wframe = None
# wframe = ax.plot_wireframe(X, Y, groups[390], rstride = 1, cstride = 1)
# print(groups[390])
# plt.show()

for i in range(0,time):
	oldcol = wframe

	Z = groups[i]
	wframe = ax.plot_wireframe(X, Y, Z, rstride = 1, cstride = 1)

	if oldcol is not None:
		ax.collections.remove(oldcol)

	plt.pause(.1)

plt.show()