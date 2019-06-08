import numpy as np
import matplotlib.pylab as pl

x = np.loadtxt('./data/x.txt')
y = np.loadtxt('./data/y.txt')
p = np.loadtxt('./data/p.txt')
print(sum(p))
pl.figure(figsize=(10.,10.), dpi=1000)
pl.xlim((300,500))
pl.ylim((250,450))
pl.tricontourf(x,y,p, 300,cmap='hot')
#pl.tricontour(x,y,p, 300,cmap='hot')
#pl.scatter(x,y,c=p, cmap = 'hot')
pl.colorbar()
pl.savefig('TW.png')