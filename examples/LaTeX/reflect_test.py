from sympy import *
from ga import Ga
from mv import MV
from printer import Format, xpdf, Fmt
Format()

ew,ex,ey,ez = MV.setup('e_w e_x e_y e_z',metric=[1,1,1,1])

a = MV('a','vector')
a.set_coef(1,0,0)
b = ex+ey+ez
c = MV('c','vector')
print 'a =',a
print 'b =',b
print 'c =',c

print a.reflect_in_blade(ex^ey).Fmt(1,'a\\mbox{ reflect in }xy')
print a.reflect_in_blade(ey^ez).Fmt(1,'a\\mbox{ reflect in }yz')
print a.reflect_in_blade(ez^ex).Fmt(1,'a\\mbox{ reflect in }zx')
print a.reflect_in_blade(ez^(ex+ey)).Fmt(1,'a\\mbox{ reflect in plane }(x=y)')
print b.reflect_in_blade((ez-ey)^(ex-ey)).Fmt(1,'b\\mbox{ reflect in plane }(x+y+z=0)')

print a.reflect_in_blade(ex).Fmt(1,'\\mbox{Reflect in }\\bm{e}_{x}')
print a.reflect_in_blade(ey).Fmt(1,'\\mbox{Reflect in }\\bm{e}_{y}')
print a.reflect_in_blade(ez).Fmt(1,'\\mbox{Reflect in }\\bm{e}_{z}')

print c.reflect_in_blade(ex^ey).Fmt(1,'c\\mbox{ reflect in }xy')
print c.reflect_in_blade(ex^ey^ez).Fmt(1,'c\\mbox{ reflect in }xyz')
print (ew^ex).reflect_in_blade(ey^ez).Fmt(1,'wx\\mbox{ reflect in }yz')
print (ew^ex).reflect_in_blade(ex^ey).Fmt(1,'wx\\mbox{ reflect in }xy')


xpdf()
