
import ga,sympy
coords=sympy.symbols('x,y,z',real=True)
base=ga.Ga('e0 e1 e2',g=[1,1,1],coords=coords)
M=[[1,2,3],[4,5,6],[7,8,9]]
A=base.lt(M)
print A
print A.lt_dict
print A.matrix()
v = base.mv('v','vector')
print v
print A(v)
