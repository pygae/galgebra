from sympy import symbols, sin
from printer import Format, xpdf, R
from ga import Ga
import sys

Format()
xyz_coords = (x, y, z) = symbols('x y z', real=True)
(o3d, ex, ey, ez) = Ga.build('e', g=[1, 1, 1], coords=xyz_coords, norm=True)
f = o3d.mv('f', 'scalar', f=True)
lap = o3d.grad*o3d.grad
#print r'%\nabla^{2} = \nabla\cdot\nabla =', lap
lap.Fmt(fmt=1,title=r'%\nabla^{2} = \nabla\cdot\nabla')
print r'%\lp\nabla^{2}\rp f =', lap*f
print r'%\nabla\cdot\lp\nabla f\rp =', o3d.grad | (o3d.grad * f)

sph_coords = (r, th, phi) = symbols('r theta phi', real=True)
(sp3d, er, eth, ephi) = Ga.build('e', g=[1, r**2, r**2 * sin(th)**2], coords=sph_coords, norm=True)
f = sp3d.mv('f', 'scalar', f=True)
lap = sp3d.grad*sp3d.grad
print r'%\nabla^{2} = \nabla\cdot\nabla =', lap
print r'%\lp\nabla^{2}\rp f =', lap*f
print r'%\nabla\cdot\lp\nabla f\rp =', sp3d.grad | (sp3d.grad * f)
print R([o3d.grad, o3d.grad])
xpdf(paper=(6, 7),debug=True)
