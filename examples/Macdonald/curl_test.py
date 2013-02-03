from sympy import symbols
from GA import *
from GAPrint import xdvi

(x,y) = symbols('x y')
(ex,ey,grad) = MV.setup('e_x e_y',metric='[1,1]',coords=(x,y))
Format()
X = MV((x,y),'vector')
print X
X3 = X*X*X
print X3
print grad*X3
print grad^(grad*X3)
xdvi()
