import sys
sys.path.append('../../')
from laga import *

basis = 'e1 e2 e3'
metric = '1 0 0, 0 1 0, 0 0 1'
(e1,e2,e3) = MV.setup(basis,metric,offset =1)

print 'e1 + (e1^e2) =',e1 + (e1^e2)
print '(e1 + (e1^e2)).is_pure() =',(e1 + (e1^e2)).is_pure()

t = symbols('t')
E = (cos(t)**2 + sin(t)**2)*e1
print E
E.trigsimp()
print E

numpy.set_printoptions(precision=3)

t,a,b = symbols('t a b')
i = e1*e2
line = cos(t)*e1 + sin(t)*e2
v = a*e1 + b*e2


A = Matrix( [ [ 1,0], [ 0,0] ] )
printeigen(A)
B = Matrix( [ [ 1,1], [ 0,0] ] )
printeigen(B)

print 'proj =', proj(line,e2)
print 'rot =', rot(i*t, e2)
print 'refl =', refl(line,e1)

M = Matrix([[1,2],[2,1]])  # Keep this
printeigen(M)

L = [Matrix(([1,2])), Matrix(([3,4]))]
GSL = GramSchmidt(L,True)
printGS(GSL)
print GramSchmidt([Matrix([1,2]), Matrix([2,1])])

u = Matrix([1,2,3])
print u.norm()

M = Matrix([ [1,1],[1,1],[0,0] ])
U, s, Vt = linalg.svd(M)
print U
print s
print Vt

A = Matrix([[ 0, 1], [ 1, 1], [ 2, 1], [3, 1]])
b = Matrix([[-1], [0.2], [0.9], [2.1]])
mhat, bhat = linalg.lstsq(A, b)[0]
print mhat, bhat

