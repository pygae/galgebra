import sys
from sympy import symbols, sin, latex, diff, Function, expand
from ga import Ga
from lt import Mlt
from printer import Eprint, Format, xpdf

Format()

#Define spherical coordinate system in 3-d

coords = (r, th, phi) = symbols('r,theta,phi', real=True)

sp3d = Ga('e_r,e_th,e_ph', g=[1, r**2, r**2*sin(th)**2], coords=coords)
(er, eth, ephi) = sp3d.mv()

#Define coordinates for 2-d (u,v) and 1-d (s) manifolds

u,v,s,alpha = symbols('u v s alpha',real=True)

sub_coords = (u,v)

smap = [1, u, v]  # Coordinate map for sphere of r = 1 in 3-d

print(r'(u,v)\rightarrow (r,\theta,\phi) = ',smap)

#Define unit sphere manifold

sph2d = sp3d.sm(smap,sub_coords)

print('#Unit Sphere Manifold:')

print('g =',sph2d.g)

(eu,ev) = sph2d.mv()

#Define vector and vector field on unit sphere tangent space

a = sph2d.mv('a','vector')
b = sph2d.mv('b','vector')
c = sph2d.mv('c','vector')
f = sph2d.mv('f','vector',f=True)

print('a =', a)
print('f =', f)
print(r'%\nabla =', sph2d.grad)

#Define directional derivative in direction a for unit sphere manifold

dd = a|sph2d.grad

print(r'%a\cdot\nabla =', dd)
print(r'%\paren{a\cdot\nabla}\bm{e}_u =', dd * eu)
print(r'%\paren{a\cdot\nabla}\bm{e}_v =', dd * ev)
print(r'%\paren{a\cdot\nabla}f =',dd * f)
print(r'%\nabla f =', sph2d.grad * f)

a1 = sph2d.mv('a_1','vector')
a2 = sph2d.mv('a_2','vector')

#Define curve on unit sphere manifold

us = Function('u__s')(s)
vs = Function('v__s')(s)

#Define 1-d submanifold on unit shpere manifold

crv1d = sph2d.sm([us,vs],[s])

(es,) = crv1d.mv()

print('#1-D Manifold On Unit Sphere:')

print(r'%\nabla =', crv1d.grad)

#Define scalar and vector fields on 1-d manifold tangent space

g = crv1d.mv('g','scalar',f=True)
h = crv1d.mv('h','vector',f=True)

print(r'%\nabla g =', crv1d.grad * g)
print(r'%\nabla \cdot \bm{h} =', crv1d.grad | h)

xpdf()
