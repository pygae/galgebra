from __future__ import print_function

import sympy
from galgebra import ga

coords=sympy.symbols('x,y,z',real=True)
base=ga.Ga('e0 e1 e2',g=[1,1,1],coords=coords)
M=[[1,2,3],[4,5,6],[7,8,9]]
A=base.lt(M)
print(A)
e0,e1,e2=base.basis
print('A.lt_dict[e0]=', A.lt_dict[e0])
print('A.lt_dict[e1]=', A.lt_dict[e1])
print('A.lt_dict[e2]=', A.lt_dict[e2])
print(A.matrix())
v = base.mv('v','vector')
print(v)
print(A(v))
