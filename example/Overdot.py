
from sympy import symbols, sin, cos
from ga import Ga
from mv import Mv, Pdop, Sdop
from printer import Format, xpdf
Format()

coords = (x,y) = symbols('x y',real=True)

o2d = Ga('e_x e_y',g=[1,1],coords=coords)
grad2d = o2d.grad
A = o2d.mv('A','mv',f=True)
B = o2d.mv('B','mv',f=True)
print(r'\bs{A}=',A)
print(r'\bs{B}=',B)
print(r'\bs{AB}=',(A*B).Fmt(3))
print(r'\dot{\nabla}\bs{A} =',grad2d.odot()*A)
print(r'\bs{A}\dot{\nabla} =',A*grad2d.odot())

diff = grad2d*(A*B) - (grad2d*A)*B -(grad2d.odot()*A)*B.odot()
print(r'\nabla\lp\bs{AB}\rp-\lp\nabla\bs{A}\rp\bs{B}-\lp\dot{\nabla}\bs{A}\rp\dot{\bs{B}}=',diff)
diff = (A*B).odot()*grad2d.odot() - A*(B.odot()*grad2d.odot()) - A.odot()*(B*grad2d.odot())
print(r'\lp\bs{\dot{AB}}\rp\dot{\nabla}-\bs{A}\lp\dot{\bs{B}}\dot{\nabla}\rp-\dot{\bs{A}}\lp\bs{B}\dot{\nabla}\rp=',diff)
C = A^B

print(r'\bs{A\W B} =',C)

diff = (grad2d^C) - ((grad2d^A)^B) - ((grad2d.odot()^A)^B.odot())
print(r'\nabla\W\lp\bs{A\W B}\rp - \lp\nabla\W \bs{A}\rp\W \bs{B} - \lp\dot{\nabla}\W A\rp\W\dot{\bs{B}}=',diff)

diff = (C.odot()^grad2d.odot()) - (A^(B.odot()^grad2d.odot())) - (A.odot()^(B^grad2d.odot()))
print(r'\lp\dot{\bs{A\W B}}\rp\W\dot{\nabla} - \bs{A} \W\lp\dot{\bs{B}}\W\dot{\nabla}\rp - \dot{\bs{A}}\W\lp\bs{B}\W\dot{\nabla}\rp=',diff)

xpdf(paper=(28,12))
