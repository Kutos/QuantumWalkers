import numpy as np
import matplotlib.pylab as pl

x = np.linspace(-100,99,200)
p1 = np.loadtxt('p.txt')

pl.figure(0,figsize=(5.,5.), dpi=300)
pl.plot(x,p1,label='1/sqrt(2) * (|up>+i|down>) at t=0')
pl.legend()
pl.savefig('HW.pdf')