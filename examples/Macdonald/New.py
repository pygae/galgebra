
from sympy import sin,cos
from GA import *

def expbi(itheta):    # Expand e^(itheta) = cos(theta) + isin(theta)
    theta = itheta.norm()
    i = itheta/theta
    print i
    result = cos(theta) + i*sin(theta)
    #result.trigsimp()  # trigsimp not working well
    return( result )

basis  = 'e_1 e_2 e_3'
metric = '1 0 0, 0 1 0, 0 0 1'

(e_1,e_2,e_3) = MV.setup(basis,metric,rframe=True)

X = e_1^e_2
print X
Y = e_1+e_3
print Y
print rot(X,Y).simplify()









