from sympy import symbols, sin
from printer import Format, xpdf, Fmt
from ga import Ga

Format()

g = '# 0 #, 0 # 0, # 0 #'

(g3d, ea, eab, eb) = Ga.build('e_a e_ab e_b', g=g)

print(g3d.g)

v = g3d.mv('v','vector')
B = g3d.mv('B','bivector')

print(v)
print(B)
xpdf()
