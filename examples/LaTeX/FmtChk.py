from sympy import symbols, sin
from printer import Format, xpdf, Fmt
from ga import Ga
import sys

Format()
xyz_coords = (x, y, z) = symbols('x y z', real=True)
(o3d, ex, ey, ez) = Ga.build('e', g=[1, 1, 1], coords=xyz_coords, norm=True)
f = o3d.mv('f', 'scalar', f=True)
F = o3d.mv('F', 'vector', f=True)
B = o3d.mv('B', 'bivector', f=True)
l = [f,F,B]
print Fmt(l)
print Fmt(l,1)
print F.Fmt(3)
print B.Fmt(3)

lap = o3d.grad*o3d.grad
print r'%\nabla^{2} = \nabla\cdot\nabla =', lap
dop = lap + o3d.grad
print dop.Fmt(fmt=3,dop_fmt=3)

xpdf(paper=(6, 7))
