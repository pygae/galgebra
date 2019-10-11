from galgebra.ga import Ga
from galgebra.printer import Format, xtex

def cross(x, y):
    return (x ^ y).dual()

Format()
GA = Ga('e*1|2|3')
print('g =',GA.g)
a = GA.mv('a', 'vector')
b = GA.mv('b', 'vector')
c = GA.mv('c', 'vector')
(e1,e2,e3) = GA.mv()
print('a =',a)
print('b =',b)
print('c =',c)

print('I =',GA.i)
e123 = e1^e2^e3
print(r'(e_\W e_2\W e_3)^2 =',e123*e123)

bc = e123*(b^c)
print(r'(e_\W e_2\W e_3)(b\W c) =',bc.Fmt(3))
abc = a^bc
print(r'a\W ((e_\W e_2\W e_3)(b\W c)) =',abc.Fmt(3))

xtex(paper=(40,11))
