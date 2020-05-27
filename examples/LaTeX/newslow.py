from galgebra.ga import Ga
from sympy import Symbol
from galgebra.printer import Format, xtex
Format()
GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

e_0_inv = e_0.inv()

print('e_0^{-1} =', e_0_inv)


p = GA.mv((1, Symbol('p1'), Symbol('p2'), Symbol('p3')), 'vector')
q = GA.mv((0, Symbol('q1'), Symbol('q2'), Symbol('q3')), 'vector')
r = GA.mv((0, Symbol('r1'), Symbol('r2'), Symbol('r3')), 'vector')
print('p =',p)
print('q =',q)
print('r =',r)

A = q ^ r
X = (p ^ A)                                 # A plane

print(r'q\W r = A =',A)
print(r'p\W A = X =',X)

D =  e_0_inv<X

print(r'e_0^{-1}\lfloor X = D =', e_0_inv<X)
#Dinv = D.rev()/(D*D.rev()).scalar()
Dinv = D.inv()
print('D^{-1} =', Dinv.Fmt(3))
print('DD^{-1} =', D*Dinv)
N = e_0_inv < (e_0 ^ X)
print(r'e_0^{-1}\lfloor (e_0\W X) = N', N)

d = (e_0_inv < (e_0 ^ X)) * Dinv  # Support vector from plane

print(r'(e_0^{-1}\lfloor (e_0\W X))  D^{-1} =',d.Fmt(3))

print('N/D =',N/D)
xtex()
