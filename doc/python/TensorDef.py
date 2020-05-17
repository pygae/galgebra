
import sys
from sympy import symbols, sin, cos
from galgebra.printer import Format, xpdf, Print_Function
from galgebra.ga import Ga
from galgebra.lt import Mlt

coords = symbols('t x y z', real=True)
(st4d, g0, g1, g2, g3) = Ga.build('gamma*t|x|y|z', g=[1, -1, -1, -1],
                                  coords=coords)

A = st4d.mv('T', 'bivector')


def TA(a1, a2):
    return A | (a1 ^ a2)


T = Mlt(TA, st4d)  # Define multi-linear function
