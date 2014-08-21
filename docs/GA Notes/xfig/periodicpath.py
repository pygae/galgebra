#!/usr/bin/env python
"""
Example: simple line plot.
Show how to make a plot that has equal aspect ratio
"""
from matplotlib import rc
from pylab import *

def rho(tau,tau0):
    sintau = sin(tau)
    csctau = 1.0/sintau
    costau = cos(tau)
    cottau = costau/sintau
    sintau0 = sin(tau0)
    csctau0 = 1.0/sintau0
    costau0 = cos(tau0)
    cottau0 = costau0/sintau0
    return(sintau*log(abs((csctau0+cottau0)/(csctau+cottau))))

rc('text', usetex=True)

tau0 = 0.00001

tau = arange(tau0,pi-tau0,0.01)

plot(tau,rho(tau,tau0),lw=2)

xlabel(r'$\frac{\tau}{T}$',fontsize=24)
ylabel(r'$\frac{\rho\left ( \tau\right )}{T}$',fontsize=24)
title(r'Photon Path $\left ( r\left ( \tau\right) = \sin\left ( \frac{\tau}{T}\right)\right)$')
#legend(loc='upper left')


grid(True)

#axes().set_aspect('equal', 'datalim')


show()
