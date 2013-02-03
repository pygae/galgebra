from sympy import *


M = Matrix([ [1, .6, .6], [.6, 1, .9], [.6, .9, 1] ])
#M = Matrix([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
M = M.evalf(6)
print 'M =', M
print 'eig values =',M.eigenvals()
print 'eig vectors =',M.eigenvects()


