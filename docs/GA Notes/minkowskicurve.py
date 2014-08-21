#!/usr/bin/env python
"""
Example: simple line plot.
Show how to make a plot that has equal aspect ratio
"""
from matplotlib import rc
from matplotlib.pyplot import quiver,quiverkey
import matplotlib.pyplot as plt
import numpy as np
import math
from pylab import *

def timeline():
    def f(tau):
        return(tau**2/2.0)
    def g(tau):
        return(0.5*(tau*math.sqrt(1+tau**2)+math.asinh(tau)))
    def df(tau):
        return(tau)
    def dg(tau):
        return(math.sqrt(1+tau**2))
    return(f,g,df,dg)

rc('text', usetex=True)
(f,g,df,dg) = timeline()
x = np.vectorize(f)
y = np.vectorize(g)
dx = np.vectorize(df)
dy = np.vectorize(dg)
tau = np.arange(0.0,4,0.01)

xlabel(r'$r( \tau )$',fontsize=24)
ylabel(r'$t( \tau )$',fontsize=24)
#legend(loc='upper left')
plot(x(tau),y(tau),lw=2)
grid(True)
TAU = np.arange(0.01,4,0.2)
lpoints = plot(x(TAU),y(TAU), 'ro')
setp(lpoints, 'markersize', 8)
setp(lpoints, 'markerfacecolor', 'r')
axes().set_aspect('equal', 'datalim')


show()
