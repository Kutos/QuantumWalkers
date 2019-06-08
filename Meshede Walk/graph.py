# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pl
import numpy as np

x = np.loadtxt('x.txt')
y = np.loadtxt('y.txt')
p = np.loadtxt('p.txt')

fig1 = pl.figure(1);
ax1 = fig1.gca(projection='3d')
ax1.plot_surface(x, y, p, cmap = 'viridis', linewidth=0.2, antialiased=True)
ax1.legend()