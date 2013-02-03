from sympy import *
from GA import *
from Manifold import *

enhance_print()

basis = 'e_x e_y e_z'
metric = '1 0 0, 0 1 0, 0 0 1'
coords = (x,y,z) = symbols('x y z')
(ex, ey, ez, grad) = MV.setup(basis,metric,coords)

f = x**2*cos(y)*ex + y*log(z)*ey + (y+1)*log(x)*ez
f = x**2 + 2*y**2

(h1, h2, h3) = symbols('h_1 h_2 h_3')
h = h1*ex + h2*ey + h3*ez

print 'Directional derivative =', DD(h,f)

# Define surface
mfvar = (u,v) = symbols('u v')
X = u*ex+v*ey+(u**2+v**2)*ez
MF = Manifold(X,mfvar)
(eu,ev) = MF.basis
# Define field on the surface.
g = (v+1)*log(u)
print 'g =',g
print 'MF.DD(eu+ev,g) =',MF.DD(eu+ev,g,True)
print 'MF.DD(eu+ev,g) =',MF.DD(eu+ev,g)
