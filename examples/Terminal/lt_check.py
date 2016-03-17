import ga,sympy
base=ga.Ga('e0 e1',g=[1,1],coords=sympy.symbols('x,y',real=True))
e0,e1 = base.mv()
M=[[1,2],[3,4]]
Mvec = [e0+2*e1,3*e0+4*e1]
#print M[0][0],M[0][1]
#print M[1][0],M[1][1]
print base.lt(Mvec)
print base.lt(M)
print sympy.Matrix(M)
