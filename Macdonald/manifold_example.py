from sympy import *
from mv import Mv
from ga import Ga
from printer import Format, xpdf

Format()

u,v = symbols('u v',real=True)

M,eu,ev = Ga.build('e_u e_v',coords=(u,v),X=(u,v,u**2+v**2))

grad = M.grad
f = (v+1)*eu + u**2*ev
print r'\Nabla =',grad
print 'g =',M.g

Du = (eu|grad).simplify()
Dv = (ev|grad).simplify()

print r'\bm{e}_{u}\cdot grad =',Du
print r'\bm{e}_{v}\cdot grad =',Dv

print r'\partial_{u}\bm{e}_{u} =',Du*eu
print r'\partial_{v}\bm{e}_{u} =',Dv*eu
print r'\partial_{u}\bm{e}_{v} =',Du*ev
print r'\partial_{v}\bm{e}_{v} =',Dv*ev

print 'f =',f

print 'grad * f =',grad*f
print 'grad | f =',grad|f
print 'grad ^ f =',grad^f

h = eu + ev
hd = h|grad

print r'((\eb_{u}+\eb_{v}) \cdot grad)f =',hd*f

r,th = symbols('r theta',real=True)

Mrth,er,eth = Ga.build('e_r e_theta',coords=(r,th),X=(r*cos(th),r*sin(th),r**2))

grad_rth = Mrth.grad
f = (r+1)*er + r**2*eth
print r'\Nabla =',grad_rth
print 'g =',Mrth.g

Du = (er|grad_rth).simplify()
Dv = (eth|grad_rth).simplify()

print r'\bm{e}_{u}\cdot grad =',Du
print r'\bm{e}_{v}\cdot grad =',Dv

print r'\partial_{u}\bm{e}_{u} =',Du*er
print r'\partial_{v}\bm{e}_{u} =',Dv*er
print r'\partial_{u}\bm{e}_{v} =',Du*eth
print r'\partial_{v}\bm{e}_{v} =',Dv*eth

print 'f =',f

print 'grad * f =',grad_rth*f
print 'grad | f =',grad_rth|f
print 'grad ^ f =',grad_rth^f

h = er + eth
hd = h|grad_rth

print r'((\eb_{r}+\eb_{\theta}) \cdot grad)f =',hd*f

xpdf(paper='landscape')

