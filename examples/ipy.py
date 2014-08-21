from printer import  Format, xpdf, Eprint
from ga import Ga
from sympy import symbols

Format()
#Eprint()

coords = (x,y,z) = symbols('x,y,z',real=True)
(o3d,ex,ey,ez) = Ga.build('e_x e_y e_z',g=[1,1,1],coords=coords)

v = o3d.mv('v','vector')
v.Fmt(3,'v')
V = o3d.mv('V','vector',f=True)

V.Fmt()

gradV = o3d.grad*V
gradV.Fmt(3,r'\nabla V')
gradV.Fmt(2,r'\nabla V')

lap = o3d.grad|o3d.grad
lap.Fmt(1,'%\\nabla^{2}')
xpdf()
