import numpy as np
import matplotlib.pylab as pl

x = np.loadtxt('./data/x.txt')
y = np.loadtxt('./data/y.txt')
p = np.loadtxt('./data/p1.txt')
print(sum(p))
pl.figure()
pl.xlim((800,1300))
pl.ylim((750,1150))
pl.tricontourf(x,y,p, 300,cmap='hot')
pl.tricontour(x,y,p, 300,cmap='hot')
pl.colorbar()
pl.show()

p = np.loadtxt('./data/p2.txt')
print(sum(p))
pl.figure()
pl.xlim((800,1300))
pl.ylim((750,1150))
pl.tricontourf(x,y,p, 300,cmap='hot')
pl.tricontour(x,y,p, 300,cmap='hot')
pl.colorbar()
pl.show()