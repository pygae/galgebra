 from sympy import *
from GA import *
from GAprecedence import GAeval,define_precedence

define_precedence(globals(),'<>|,*,^') #Macdonald precedence

basis = 'ex ey ez'
metric = '1 0 0, 0 1 0, 0 0 1'
coords = (x, y, z) = symbols('x y z')
(ex,ey,ez) = MV.setup('ex ey ez',metric)

print GAeval('ex + 4*ex ^ ey < ex*ez',True)



