#!/usr/bin/python

from sympy import *
from GA import *
from GAPrint import xdvi

Format()

rho, phi = symbols('rho phi')
X = Matrix([ [rho*cos(phi)], [rho*sin(phi)] ])
Y = Matrix([ [rho, phi] ])
#X is a 1 x m matrix of functions. Y is a 1 x n matrix of variables to differentiate with respect to.
#Then X.jacobian(Y) is a m x n matrix, the differential of X.

print 'X =',X
print 'Y =',Y

Diff = X.jacobian(Y)
print '\\pdiff{X}{Y} =',Diff

print '\\abs{\\pdiff{X}{Y}} = ', trigsimp(Diff.det())

DiffInv = Diff.inv()
DiffInv.simplify()
print '%\\lp\\pdiff{X}{Y}\\rp^{-1} =',DiffInv

xdvi()

