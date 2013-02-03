from GA import *
from GAPrint import enhance_print,get_program,print_function
from sympy import *
from Manifold import *

enhance_print()

#Format('1 1 1 1')
Format()

#Format(1)
coords = (x,y,z) = symbols('x y z')
basis = (ex, ey, ez, grad) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=coords)


# Define surface
mfvar = (u,v) = symbols('u v')
X = u*ex+v*ey+(u**2+v**2)*ez
X.Fmt(1,'Manifold')
MF = Manifold(X,mfvar)
(eu,ev) = MF.Basis()

# Define field on the surface.
g = (v+1)*log(u)

dg = MF.grad*g
dg.Fmt(3,'dg')
dg.subs({u:1,v:0}).Fmt(1,'dg(1,0)')

# Define vector field on the surface

G = v**2*eu+u**2*ev
G.Fmt(1,'G')
dG = MF.grad*G
dG.Fmt(3,'dG')
dG.subs({u:1,v:0}).Fmt(1,'dG(1,0)')

xdvi()
