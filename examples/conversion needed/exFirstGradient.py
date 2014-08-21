from sympy import *
from ga import *
from ga_print import xdvi

Format()
basis = 'e1 e2 e3'
metric = '[1 ,1 ,1]'
coords = (x,y,z) = symbols('x y z')
(e1,e2,e3) = MV.setup('e_1 e_2 e_3', metric)
(e1, e2, e3, grad) = MV.setup(basis,metric,coords)

# Define f and take its gradient
f = 1 + x*y^2*e1 + sin(y)*e2 + sin(x)*cos(y)*e1*e2

print 'f =', f
print 'grad f =', grad * f

f = 1 + (x*y**2)*e1 + sin(y)*e2 + sin(x)*cos(y)*e1*e2

print 'f =', f
print 'grad f =', grad * f

print r'grad \cdot f =', grad < f
print 'grad \W f =', grad ^ f

xdvi()
