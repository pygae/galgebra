from __future__ import print_function
import sys
from galgebra.printer import  Format, xpdf, xtex
Format()
from sympy import symbols, sin, pi, latex, Array, permutedims
from galgebra.ga import Ga

# From http://latexcolor.com/
print(r'\definecolor{airforceblue}{rgb}{0.36, 0.54, 0.66}')
print(r'\definecolor{applegreen}{rgb}{0.55, 0.71, 0.0}')
print(r'\definecolor{atomictangerine}{rgb}{1.0, 0.6, 0.4}')

print(r'\bm{\mbox{Base manifold (three dimensional)}}')
print(r'\bm{\mbox{Metric tensor (cartesian coordinates - norm = False)}}')
from sympy import cos, sin, symbols
g3coords = (x,y,z) = symbols('x y z')
g3 = Ga('ex ey ez', g = [1,1,1], coords = g3coords,norm=False) # Create g3
(e_x,e_y,e_z) = g3.mv()

print('g =',g3.g)
print('\\')

print(r'\bm{\mbox{Two dimensioanal submanifold - Unit sphere}}')
print(r'\text{Basis not normalised}')

sp2coords = (theta, phi) = symbols(r'theta phi', real = True)
sp2param = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]

sp2 = g3.sm(sp2param, sp2coords, norm = False) # submanifold

(etheta, ephi) = sp2.mv() # sp2 basis vectors
(rtheta, rphi) = sp2.mvr() # sp2 reciprocal basis vectors

sp2grad = sp2.grad

sph_map = [1, theta, phi]  # Coordinate map for sphere of r = 1
print(r'(\theta,\phi)\rightarrow (r,\theta,\phi) = ', latex(sph_map))

(etheta, ephi) = sp2.mv()
print(r'e_\theta | e_\theta = ', etheta|etheta)
print(r'e_\phi | e_\phi = ', ephi|ephi)

print('g =',sp2.g)
print(r'\text{g\_inv = }', latex(sp2.g_inv))

#print(r'\text{signature = ', latex(sp2.signature()))

Cf1 = sp2.Christoffel_symbols(mode=1)
Cf1 = permutedims(Array(Cf1), (2, 0, 1))
print(r'\text{Christoffel symbols of the first kind: }')
print(r'\Gamma_{1, \alpha, \beta} = ', latex(Cf1[0, :, :]), r'\quad', r'\Gamma_{2, \alpha, \beta} = ', latex(Cf1[1, :, :]))

Cf2 = sp2.Christoffel_symbols(mode=2)
Cf2 = permutedims(Array(Cf2), (2, 0, 1))
print(r'\text{Christoffel symbols of the second kind: }')
print(r'\Gamma^{1}_{\phantom{1,}\alpha, \beta} = ', latex(Cf2[0, :, :]), r'\quad', r'\Gamma^{2}_{\phantom{2,}\alpha, \beta} = ', latex(Cf2[1, :, :]))

F = sp2.mv('F','vector',f=True) #scalar function)
f = sp2.mv('f','scalar',f=True) #vector function)
print(r'\nabla = ', sp2grad)
print(r'\nabla f =',sp2.grad * f)
print('F =',F)
print(r'\nabla F = ',sp2.grad * F)
print('\\')

print(r'\mbox{One dimensioanal submanifold}')
print(r'\mbox{Basis not normalised}')

#cir_th = phi = symbols(r'{\color{atomictangerine}\phi}',real = True)
cir_th = phi = symbols('phi',real = True)
cir_map = [pi/8, phi]
print(r'(\phi)\rightarrow (\theta,\phi) = ', latex(cir_map))

cir1d = sp2.sm( cir_map , (cir_th,), norm = False) # submanifold

cir1dgrad = cir1d.grad

(ephi) = cir1d.mv()
print(r'e_\phi | e_\phi = ', latex(ephi[0] | ephi[0]))
print('g = ', latex(cir1d.g))

h = cir1d.mv('h','scalar',f= True)

H = cir1d.mv('H','vector',f= True)

print(r'\nabla = ', cir1dgrad)
print(r'\nabla h = ', (cir1d.grad * h).simplify())
print('H =', H)
print(r'\nabla H = ', (cir1d.grad * H).simplify())
print('\\' )

# xpdf(paper=(9,10))
#xpdf(paper=(9,10),pdfprog=None)
xtex()

