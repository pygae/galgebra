# Provided by Alan Macdonald   http://faculty.luther.edu  macdonal@luther.edu

from sympy import *
from ga import *

# g2: Standard 2D Model
g2coords = (x,y) = symbols('x y', real=True)
g2 = Ga('e_x e_y', g=[1,1], coords=g2coords)
(ex, ey) = g2.mv()

# g3: Standard 3D Model
g3coords = (x,y,z) = symbols('x y z', real=True)
g3 = Ga('e_x e_y e_z', g=[1,1,1], coords=g3coords)
(ex, ey, ez) = g3.mv()

# g4: Standard 4D Model
g4coords = (w,x,y,z) = symbols('w x y z', real=True)
g4 = Ga('e_w e_x e_y e_z', g=[1,1,1,1], coords=g4coords)
(ew, ex, ey, ez) = g4.mv()

# h3: Homogeneous Coordinates
h3coords = (x,y,z,e) = symbols('x y z e', real=True)
h3 = Ga('e_x e_y e_z e_e', g=[1,1,1,1], coords=h3coords)
(ex,ey,ez,ee) = h3.mv()

# sp3: Spherical Coordinates I
sp3coords = (r, ph, th) = symbols('r phi theta', real=True)
sp3 = Ga('e_r e_phi e_theta', g=None, coords=sp3coords, \
          X=[r, r*sin(ph)*cos(th), r*sin(ph)*sin(th), r*cos(ph)], norm=True)
(er, eph, eth) = sp3.mv()
# Mathematics naming convention: $\phi$ colatitude, $\theta$ longitude.
# (r phi theta) is a right-handed system

# sp3: Spherical Coordinates II
sp3coords = (r, ph, th) = symbols('r phi theta', real=True, norm=True)
sp3 = Ga('e_r e_phi e_theta', g=[1, r**2, r**2*sin(ph)**2], coords=sp3coords)
(er, eph, eth) = sp3.mv()
# Mathematics naming convention: $\phi$ colatitude, $\theta$ longitude.
# (r phi theta) is a right-handed system

# st4: Spacetime Algebra I
st4coords = (t,x,y,z) = symbols('t x y z', real=True)
st4 = Ga('e_t e_x e_y e_z', g=[1,-1,-1,-1], coords=st4coords)
(et, ex, ey, ez) = st4.mv()

# st4: Spacetime Algebra II
st4coords = (t,x,y,z) = symbols('t x y z', real=True)
st4 = Ga('e_t e_x e_y e_z', g=[-1,1,1,1], coords=st4coords)
(et, ex, ey, ez) = st4.mv()

# cf3: Conformal Model I  (Doran/Lasenby p. 363)
cf3coords = (n,x,y,z,nbar) = symbols('n x y z nbar', real=True)
cf3g = '0 0 0 0 2, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, 2 0 0 0 0'
cf3 = Ga('e_n e_x e_y e_z e_nbar', g = cf3g, coords = cf3coords)
(en, ex, ey, ez, enbar) = cf3.mv()

# sp2: Unit Sphere in $\RE^3$ I  (submanifold of g3)
sp2coords = (ph,th) = symbols('phi theta', real=True)
sp2param = [sin(ph)*cos(th), sin(ph)*sin(th), cos(ph)]
sp2 = g3.sm(sp2param, sp2coords, norm=True)
(eph, eth) = sp2.mv()

# sp2: Unit Sphere in $\RE^3$ II  (submanifold of sp3)
sp2coords = (p,t) = symbols('phi theta', real=True)  # p for phi, t for theta
sp2param = [1, p, t]  # Parameterization of unit sphere
sp2 = sp3.sm(sp2param, sp2coords, norm=True)
(ep, et) = sp2.mv()
