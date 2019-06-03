import numpy as np
import matplotlib.pylab as pl

x = np.loadtxt('./data/x.txt')
y = np.loadtxt('./data/y.txt')
p = np.loadtxt('./data/p.txt')
print(sum(p))
pl.tricontourf(x,y,p, 300,cmap='hot')
pl.tricontour(x,y,p, 300,cmap='hot')
pl.savefig('nstep500.pdf')