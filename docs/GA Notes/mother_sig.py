from sympy import *

A = Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])

print A

print A.eigenvals()

A = Matrix([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])

print A

print A.eigenvals()

A = Matrix([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[1,0,0,0,-1]])

print A

print A.eigenvals()

A = Matrix([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,1],[1,0,0,1,0]])

print A

print A.eigenvals()
