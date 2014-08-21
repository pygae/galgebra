#!/usr/bin/env python
"""
Example: simple line plot.
Show how to make a plot that has equal aspect ratio
"""
from matplotlib import rc
from pylab import *

def rho(tau,eta):
    if eta == 1.0:
        return(tau*log(tau))
    else:
        return((1.0/(1.0-float(eta))*(tau-tau**float(eta))))

rc('text', usetex=True)

tau = arange(1.0,3.0,0.01)
plot(tau,rho(tau,-1.5),lw=2,label='-1.5')
plot(tau,rho(tau,-1.0),lw=2,label='-1.0')
plot(tau,rho(tau,-0.5),lw=2,label='-0.5')
plot(tau,rho(tau,0.0),lw=2,label=' 0.0')
plot(tau,rho(tau,0.5),lw=2,label=' 0.5')
plot(tau,rho(tau,1.0),lw=2,label=' 1.0')
plot(tau,rho(tau,1.5),lw=2,label=' 1.5')
xlabel(r'$\frac{\tau}{\tau_{0}}$',fontsize=24)
ylabel(r'$\frac{\rho\left ( \tau\right )}{\tau_{0}}$',fontsize=24)
title(r'Photon Path $\left ( r\left (\tau \right ) = \left ( \frac{\tau}{\tau_{0}}\right )^{\eta}\right )$')
legend(loc='upper left')


grid(True)

#axes().set_aspect('equal', 'datalim')


show()
