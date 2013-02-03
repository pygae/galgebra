from GAPrint import xdvi
from GA import *

Format()
X = (x,y,z) = symbols('x y z')
(ex,ey,ez,grad) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=X)

f = MV('f','scalar',fct=True)
A = MV('A','vector',fct=True)
B = MV('B','grade2',fct=True)

print r'\bm{A} =',A
print r'\bm{B} =',B

print 'grad*f =',grad*f
print r'grad|\bm{A} =',grad|A
(grad*A).Fmt(2,r'grad*\bm{A}')

print r'-I*(grad^\bm{A}) =',-MV.I*(grad^A)
(grad*B).Fmt(2,r'grad*\bm{B}')
print r'grad^\bm{B} =',grad^B
print r'grad|\bm{B} =',grad|B

xdvi(debug=True)
