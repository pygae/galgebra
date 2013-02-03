from sympy import *
from GA import *
#Python precedence: *, +-, ^, |, <>             # GOI
#define_precedence(globals(), op_ord='<>|,^,*') # IOG  # Default is Doran and Lasenby convention
define_precedence(globals(), '*,<>|,^')         # GIO  # Macdonald convention.
from Manifold import *

#enhance_print()
Format()

coords = (x,y,z) = symbols('x y z')
(ex,ey,ez,grad) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=coords)

# Define a manifold M
Mvar = (u,v) = symbols('u v')
X = u*ex + v*ey + u*v*ez
M = Manifold(X,Mvar)

# Define a field f on M.
(xu,xv) = M.basis
f = MV(u**2+v**2,'scalar')

print 'grad*f ^ V = A|B'

xdvi()
