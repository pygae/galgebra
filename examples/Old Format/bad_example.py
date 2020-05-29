from __future__ import print_function
from sympy import S
from galgebra.printer import xdvi,Print_Function,Format
from galgebra.deprecated import MV

(ex,ey,ez) = MV.setup('e_x e_y e_z',metric='[1,1,1]')

v = S(3)*ex+S(4)*ey

print(v)
print(v.norm())
print(v.norm2())
print(v/v.norm())
print(v/(v.norm()**2))
print(v*(v/v.norm2()))
print(v.inv())
