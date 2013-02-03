
from sympy import symbols,sin,cos
from GAPrint import enhance_print
from GA import *

enhance_print()

X = (x,y,z) = symbols('x y z')
(ex,ey,ez,grad) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=X)

f = MV(x**2+y**2+z**2,'scalar')

print 'f =',f
print 'grad*f =',grad*f
print 'grad*f(1,2,3) =',(grad*f).subs({x:1,y:2,z:3})

A = (x**3+y**2)*ex+y**3*ey+z**3*ez

gradA = grad*A

print 'A =',A
print 'grad*A =', gradA
print 'gradA(1,2,3) =',gradA.subs({x:1,y:2,z:3})

B = z**2*(ex^ey)+y**2*(ex^ez)+x**2*(ey^ez)

print 'B =',B
print 'grad*B =',grad*B
print 'grad*B(1,2,3) =', (grad*B).subs({x:1,y:2,z:3})

print 'I^{-1} =',MV.I.inv()
