from sympy import *
from ga import *

# g3: Standard 3D Model
g3coords = (x,y,z) = symbols('x y z')
g3 = Ga('e_x e_y e_z', g=[1,1,1], coords=g3coords)
(ex, ey, ez) = g3.mv()

# g2: Standard 2D Model
g2coords = (x,y) = symbols('x y')
g2 = Ga('e_x e_y', g=[1,1], coords=g2coords)
(ex, ey) = g2.mv()

# h3: Homogeneous Coordinates
h3coords = (x,y,z,e) = symbols('x y z e')
h3 = Ga('e_x e_y e_z e_e', g=[1,1,1,1], coords=h3coords)
(ex,ey,ez,ee) = h3.mv()

# sp3: Spherical Coordinates I
sp3coords = (r, phi, th) = symbols('r  phi theta')
sp3 = Ga('e_r e_phi e_theta', g=None, coords=sp3coords, \
          X=[r, r*sin(phi)*cos(th), r*sin(phi)*sin(th), r*cos(phi)])
(er, ephi, eth) = sp3.mv()
# Mathematics naming convention: $\phi$ colatitude, $\theta$ longitude.
# (r phi theta) is a right-handed system 

# sp3: Spherical Coordinates II
sp3coords = (r, phi, th) = symbols('r phi theta', real=True)
sp3 = Ga('e_r e_phi e_theta', g=[1, r**2, r**2*sin(phi)**2], coords=sp3coords)
(er, ephi, eth) = sp3.mv()
# Mathematics naming convention: $\phi$ colatitude, $\theta$ longitude.
# (r phi theta) is a right-handed system 

# st4: Spacetime Algebra I
st4coords = (t,x,y,z) = symbols('t x y z') 
st4 = Ga('e_t e_x e_y e_z', g=[1,-1,-1,-1], coords=st4coords)
(et, ex, ey, ez) = st4.mv()

# st4: Spacetime Algebra II
st4coords = (t,x,y,z) = symbols('t x y z') 
st4 = Ga('e_t e_x e_y e_z', g=[-1,1,1,1], coords=st4coords)
(et, ex, ey, ez) = st4.mv()

# cf3g: Conformal Model I  (Doran/Lasenby p. 363)
cf3coords = (n,x,y,z,nbar) = symbols('n x y z nbar')
cf3g = '0 0 0 0 2, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, 2 0 0 0 0'
cf3 = Ga('e_n e_x e_y e_z e_nbar', g = cf3g, coords = cf3coords)
(en, ex, ey, ez, enbar) = cf3.mv()

# cfg: Conformal Model II  (Dorst p. 361)
cf3coords = (o,x,y,z,nfty) = symbols('o x y z nfty')
cf3g = '0 0 0 0 -1, 0 1 0 0 0, 0 0 1 0 0, 0 0 1 0 0, -1 0 0 0 0'
cf3 = Ga('e_o e_x e_y e_z e_nfty', g = cf3g, coords = cf3coords)
(eo, ex, ey, ez, enfty) = cf3.mv()

# sp2: Unit Sphere in $\RE^3$ I  (submanifold of g3)
sp2coords = (phi,th) = symbols('phi theta', real=True)
sp2param = [sin(phi)*cos(th), sin(phi)*sin(th), cos(phi)]
sp2 = g3.sm(sp2param, sp2coords)
(ephi, eth) = sp2.mv()

# sp2: Unit Sphere in $\RE^3$ II  (submanifold of sp3)
sp2coords = (p,t) = symbols('p t', real=True)  # p for phi, t for theta
sp2param = [1, p, t]  # Parameterization of unit sphere
sp2 = sp3.sm(sp2param, sp2coords)
(ep, et) = sp2.mv()
