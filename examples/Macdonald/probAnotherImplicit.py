from sympy import *
from GA import *

rho,theta,phi = symbols('rho theta phi')

X = Matrix([ [rho*sin(phi)*cos(theta)], [rho*sin(phi)*sin(theta)], [rho*cos(phi)] ])
Y = Matrix([ [rho, phi, theta] ])

Diff = X.jacobian(Y)
print 'Differential1 ='
print Diff

Diff_sub = Diff.subs({rho:2,phi:pi/2,theta:0})
#rho.subs(rho,2)
#phi.subs(phi,pi/2)
#theta.subs(theta,0)
#rho = 2
#phi = pi/2
#theta = 0
print "Diff_sub ="
print Diff_sub


Diff = Matrix([[1,0,0], [0,0,2], [0,-2,0]]) # So substitute by hand
print
print 'Differential ='
print Diff

print
print 'DiffInv ='
print Diff.inv()
