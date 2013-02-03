from sympy import *
from GA import *

enhance_print()

basis = 'e_x e_y e_z'
metric = '1 0 0, 0 1 0, 0 0 1'
coords = (x,y,z) = symbols('x y z')
(ex, ey, ez, grad) = MV.setup(basis,metric,coords)

F = x*ex + y*ey + (y+1)*log(x)*ez

print F.diff(x)
print F.diff(y)
print F.diff(z)
print grad*F
print grad|F
print grad^F

v = ex+2*ey+3*ez

print DD(v,F)
print MV.Iinv
