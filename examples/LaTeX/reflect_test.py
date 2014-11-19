from sympy import *
from ga import Ga
from printer import Format, xpdf, Fmt
Format()

e3d = Ga('e_x e_y e_z',g=[1,1,1])
ex,ey,ez = e3d.mv()

a = e3d.mv('a','vector')

print 'a =',a

print a.reflect_in_blade(ex^ey).Fmt(1,'\\mbox{Reflect in }xy')
print a.reflect_in_blade(ey^ez).Fmt(1,'\\mbox{Reflect in }yz')
print a.reflect_in_blade(ez^ex).Fmt(1,'\\mbox{Reflect in }zx')
print a.reflect_in_blade(ez^(ex+ey)).Fmt(1,'\\mbox{Reflect in plane }(x=y)')

print a.reflect_in_blade(ex).Fmt(1,'\\mbox{Reflect in }\\bm{e}_{x}')
print a.reflect_in_blade(ey).Fmt(1,'\\mbox{Reflect in }\\bm{e}_{y}')
print a.reflect_in_blade(ez).Fmt(1,'\\mbox{Reflect in }\\bm{e}_{z}')

xpdf()
