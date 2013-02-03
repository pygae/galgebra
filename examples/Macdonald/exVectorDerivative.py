from sympy import *
from GA import *
from Manifold import *

enhance_print()

basis = 'e_x e_y e_z'
metric = '1 0 0, 0 1 0, 0 0 1'
coords = (x,y,z) = symbols('x y z')
(ex, ey, ez, grad) = MV.setup(basis,metric,coords)

# Define surface
mfvar = (u,v) = symbols('u v')
X = u*ex + v*ey + (u**2+v**2)*ez
MF = Manifold(X,mfvar)

# Define field on the surface.
(xu,xv) = MF.basis
f = (v+1)*xu + u**2*xv

print 'Vector derivative =', MF.grad*f

