# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as pl

n = 20
g = 2 * n + 3

x = np.loadtxt('./data/x.txt')
y = np.loadtxt('./data/y.txt')
p = np.loadtxt('./data/p.txt')

print(sum(p))

p = p.reshape((g,g,g,g,3,3))

pp = np.zeros_like(x)

a = 0

for i in range(g):
    for j in range(g):
        for k1 in range(3):
            pp[a] += p[i,j,i,j,k1,k1]
            a += 1

pl.tricontourf(x,y,pp,300, cmap='hot')
#pl.scatter(x,y,c=pp)
pl.colorbar()
pl.show