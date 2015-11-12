from sympy import *
from ga import Ga

# Standard 2D Model
g2basis = 'ex ey'
g2metric = [1,1]
(x,y) = g2coords = symbols('x y',real=True)
(g2, ex, ey) = Ga.build(g2basis, g=g2metric, coords=g2coords)

# Same results with line 11 or 12.
d = g2.mv(2)
print d

print d | (3*ex)     # Crashes
print d < (3*ex)     # Output 2*ex. Should be 6*ex.
print d > (3*ex)     # Output 2*ex. Should be 0.
