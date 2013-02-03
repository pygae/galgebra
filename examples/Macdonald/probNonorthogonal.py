# Compute a reciprocal basis

from sympy import *
from GA import *

(ex,ey,ez) = MV.setup('e_x e_y e_z', metric='[1,1,1]')
#enhance_print()
Format(3)

(u,v,w) = symbols('u v w')

X = u*cos(v)*ex + u*sin(v)*ey + (w + u*cos(v))*ez
print r'\bm{x} =', X

xu = X.diff(u)
xv = X.diff(v)
xw = X.diff(w)

print r'\bm{x}_u =', xu
print r'\bm{x}_v =', xv
print r'\bm{x}_w =', xw

(xur, xvr, xwr) = ReciprocalFrame((xu,xv,xw))
print r'%\bm{x}^u =', xur
print r'%\bm{x}^v =', xvr
print r'%\bm{x}^w =', xwr


xdvi()
