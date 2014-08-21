#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)

strv = np.vectorize(str)

def F(alpha,eta):
    return(alpha**(-1.0/eta)/(np.sin(np.pi/eta)*eta))
"""
def F(alpha,eta):
    return(alpha**2-eta**2)
"""
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

delta = 0.01
eta = np.arange(1.01,1.5,delta)
alpha = np.arange(0.1,4.0,delta)
(Alpha,Eta) = np.meshgrid(alpha,eta)
Fv = np.vectorize(F)
#Fv0 = Fv(Alpha,Eta)

# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
plt.figure()
levels = np.arange(0.0,1.1,0.1)
CS = plt.contour(Alpha,Eta,Fv(Alpha,Eta),levels=levels)
plt.xlabel(r'$\alpha$',fontsize=16)
plt.ylabel(r'$\eta$',fontsize=16)
#plt.clabel(CS,inline=0, fontsize=10)
CB = plt.colorbar(CS, shrink=0.8, extend='both')
plt.title(r'$\displaystyle\frac{\alpha^{-\left (\frac{1}{\eta}\right )}}{\eta}\csc\left ( \frac{\pi}{\eta}\right )$',fontsize=16)


plt.show()
