from sympy import S
from printer import xdvi,Get_Program,Print_Function,Format
from mv import MV

(ex,ey,ez) = MV.setup('e_x e_y e_z',metric='[1,1,1]')

v = S(3)*ex+S(4)*ey

print(v)
print(v.norm())
print(v.norm2())
print(v/v.norm())
print(v/(v.norm()**2))
print(v*(v/v.norm2()))
print(v.inv())
