from __future__ import print_function
import sys

from sympy import symbols,sin,cos,latex
from galgebra.ga import Ga
from galgebra.printer import Format, xtex

g = '# 0 0 0,0 # 0 0,0 0 # 0,0 0 0 #'

Format()
(sp4d,g0,g1,g2,g3) = Ga.build('gamma*0|1|2|3',g=g)

psi = sp4d.mv('psi','even')

print(r'\psi =',psi)

B = sp4d.mv('B','bivector')

print('B =',B)

print(r'B\psi-\psi B =', (B*psi-psi*B).Fmt(3))

xtex()
