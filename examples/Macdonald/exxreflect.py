from GA import *

basis = 'e1 e2 e3'
metric = '1 0 0, 0 1 0, 0 0 1'
(e1,e2,e3) = MV.setup(basis, metric)

orig = (e1 + e3)*e2
B = e1*e2

print refl(B,orig)

print pi
