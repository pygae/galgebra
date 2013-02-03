from sympy import *
from GA import *
from GAPrint import enhance_print
enhance_print()

basis = 'e_x e_y e_z'
metric = '[1,1,1]'

MV.setup(basis,metric=metric)
#MV.setup('e*x|y|z','[1,1,1]')
Format('mv=3')

v  =  MV('V', 'vector')
print v
B  =  MV('B', 'grade2')

print B
print 'v . B =', (v*B).grade(1)
print ((v*B) + (B*v)).grade(1) # Is zero

xdvi()
