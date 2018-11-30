import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sc
from random import choice


from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 0.5 * np.outer(np.cos(u), np.sin(v))
y = 0.5 * np.outer(np.sin(u), np.sin(v))
z = 0.5 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='k', alpha=0.2)
ax.plot_surface(x*0.6/0.5, y*0.6/0.5, z*0.6/0.5, color='k', alpha=0.2)
ax.plot_surface(x*0.7/0.5, y*0.7/0.5, z*0.7/0.5, color='k', alpha=0.2)
ax.plot_surface(x*0.8/0.5, y*0.8/0.5, z*0.8/0.5, color='k', alpha=0.2)

with open("ptracks.lst") as fp:
	r = fp.read().split('&')
	hist = choice(r).split('\n')[1:-1]
	hist = [[float(x) for x in s.split(',')] for s in hist]

allpoints = [tuple(l[:3]) for l in hist]
allpoints.append(tuple(hist[-1][3:]))

allpoints = tuple(zip(*allpoints))

ax.plot(xs = allpoints[0], ys = allpoints[1], zs = allpoints[2], color='b', marker='o')

plt.show()





