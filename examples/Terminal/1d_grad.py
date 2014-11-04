from sympy import *

s = symbols('s',real=True)

x = s**2

f = Function('f')(x)

print x
print f
D_s_f = diff(f,s,evaluate=False)

print D_s_f
print latex(D_s_f)

