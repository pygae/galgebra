# Compute a reciprocal basis

from sympy import *
from GA import *
from GAPrint import enhance_print,xdvi

enhance_print()
Format(3)

coords = symbols('x y z')
(ex,ey,ez,grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]',coords=coords)

(u,v,w) = symbols('u v w')

curv = u*cos(v)*ex + u*sin(v)*ey + (w + u*cos(v))*ez

print r'\mbox{Vector Manifold}'
print r'\bm{X} =',curv

eu = curv.diff(u)
ev = curv.diff(v)
ew = curv.diff(w)

e = [eu,ev,ew]

e_str = [r'\bm{e}_{u}',r'\bm{e}_{v}',r'\bm{e}_{w}']
er_str = [r'\bm{e}^{u}',r'\bm{e}^{v}',r'\bm{e}^{w}']

print r'\mbox{Basis Vectors}'
eu.Fmt(1,e_str[0])
ev.Fmt(1,e_str[1])
ew.Fmt(1,e_str[2])

er = (eu_r, ev_r, ew_r) = ReciprocalFrame((eu,ev,ew))

print r'\mbox{Reciprocal Basis Vectors}'

eu_r.Fmt(1,'%'+er_str[0])
ev_r.Fmt(1,'%'+er_str[1])
ew_r.Fmt(1,'%'+er_str[2])

print r'\mbox{Dot Products}'

for (ei_str,ei) in zip(e_str,e):
    for (ej_str,ej) in zip(er_str,er):
        print '%'+ei_str+r'\cdot'+ej_str+' =',ei|ej

xdvi()
