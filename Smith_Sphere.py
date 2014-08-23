from ga import Ga
from printer import Format, xpdf
from sympy import *

Format()

(ga,er,ex,es) = Ga.build('e_r e_x e_s',g=[1,1,1])
(gaz,zr,zx) = Ga.build('z_r z_x',g=[1,1])

Bz = er^ex # impedance plance
Bs = es^ex # reflection coefficient plane
Bx = er^es
I = ga.I()

def down(p, N):
    '''
    stereographically project a vector in G3 downto the bivector N
    '''
    n= -1*N.dual()
    return -(n^p)*(n-n*(n|p)).inv()

def up(p):
    '''
    stereographically project a vector in G2 upto the space G3
    '''
    if (p^Bz).is_zero():
        N = Bz
    elif  (p^Bs).is_zero():
        N = Bs

    n = -N.dual()

    return   n + 2*(p*p + 1).inv()*(p-n)

a,b,c,z,s,n = [ga.mv(k,'vector') for k in ['a','b','c','z','s' ,'n']]

Bz.dual()

z = z.proj([er,ex])

p = up(z)
p.Fmt(3,'p')

p.norm2()
print down(p, Bz)
print down(p,Bs).simplify()

print (z-er)*(z+er).inv()

p.Fmt(3,'p')

R=((-pi/4)*Bx).exp()
print R

print R*p*R.rev()

print down(R*p*R.rev(),Bz)

xpdf(paper='letter')
