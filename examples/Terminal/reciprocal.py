import ga,sympy
base=ga.Ga('a b')
print(base.e)
ar,br = base.mvr()
a,b = base.mv()
print(ar)
print(br)
print((a|ar).simplify())
print(b|ar)
print(a|br)
print((b|br).simplify())
print(base.E())
print(base.I())


import ga,sympy
coords = u,v =sympy.symbols('u,v',real=True)
#x,y,z=sympy.symbols('x,y,z',Function=True)

x = sympy.Function('x')(*coords)
y = sympy.Function('y')(*coords)
z = sympy.Function('z')(*coords)

#surface=ga.Ga('e1 e2',coords=(u,v),X=[x(u,v),y(u,v),z(u,v)])
surface=ga.Ga('e1 e2',coords=(u,v),X=[x,y,z])
print(surface.g)
eu,ev = surface.mv()
eur,evr = surface.mvr()
print((eu|eur).simplify())
print((eu|evr).simplify())
print((ev|eur).simplify())
print((ev|evr).simplify())
