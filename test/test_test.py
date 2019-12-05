import sys
import unittest
from sympy import symbols, sin, cos, Rational, expand, collect, simplify, Symbol
from galgebra.printer import Format, Eprint, Get_Program, latex, GaPrinter, ZERO_STR
from galgebra.ga import Ga, one, zero
from galgebra.mv import Mv, Nga
# for backward compatibility
from galgebra.mv import ONE, ZERO, HALF
from galgebra import ga
from galgebra import utils

def F(x):
    global n, nbar
    Fx =  HALF * ((x * x) * n + 2 * x - nbar)
    return(Fx)

def make_vector(a, n=3, ga=None):
    if utils.isstr(a):
        v = zero
        for i in range(n):
            a_i = Symbol(a+str(i+1))
            v += a_i*ga.basis[i]
        v = ga.mv(v)
        return(F(v))
    else:
        return(F(a))

class TestTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_basic_multivector_operations(self):

        g3d = Ga('e*x|y|z')
        (ex,ey,ez) = g3d.mv()

        A = g3d.mv('A','mv')

        assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'
        assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'
        assert str(A) == 'A + A__x*e_x + A__y*e_y + A__z*e_z + A__xy*e_x^e_y + A__xz*e_x^e_z + A__yz*e_y^e_z + A__xyz*e_x^e_y^e_z'

        X = g3d.mv('X','vector')
        Y = g3d.mv('Y','vector')

        assert str(X) == 'X__x*e_x + X__y*e_y + X__z*e_z'
        assert str(Y) == 'Y__x*e_x + Y__y*e_y + Y__z*e_z'

        assert str((X*Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z + (X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
        assert str((X^Y)) == '(X__x*Y__y - X__y*Y__x)*e_x^e_y + (X__x*Y__z - X__z*Y__x)*e_x^e_z + (X__y*Y__z - X__z*Y__y)*e_y^e_z'
        assert str((X|Y)) == '(e_x.e_x)*X__x*Y__x + (e_x.e_y)*X__x*Y__y + (e_x.e_y)*X__y*Y__x + (e_x.e_z)*X__x*Y__z + (e_x.e_z)*X__z*Y__x + (e_y.e_y)*X__y*Y__y + (e_y.e_z)*X__y*Y__z + (e_y.e_z)*X__z*Y__y + (e_z.e_z)*X__z*Y__z'


        g2d = Ga('e*x|y')
        (ex,ey) = g2d.mv()

        X = g2d.mv('X','vector')
        A = g2d.mv('A','spinor')

        assert str(X) == 'X__x*e_x + X__y*e_y'
        assert str(A) == 'A + A__xy*e_x^e_y'

        assert str((X|A)) == '-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y)*e_x + A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y)*e_y'
        assert str((X<A)) == '-A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y)*e_x + A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y)*e_y'
        assert str((A>X)) == 'A__xy*((e_x.e_y)*X__x + (e_y.e_y)*X__y)*e_x - A__xy*((e_x.e_x)*X__x + (e_x.e_y)*X__y)*e_y'


        o2d = Ga('e*x|y', g=[1, 1])
        (ex, ey) = o2d.mv()

        X = o2d.mv('X', 'vector')
        A = o2d.mv('A', 'spinor')

        assert str(X) == 'X__x*e_x + X__y*e_y'
        assert str(A) == 'A + A__xy*e_x^e_y'

        assert str((X*A)) == '(A*X__x - A__xy*X__y)*e_x + (A*X__y + A__xy*X__x)*e_y'
        assert str((X|A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
        assert str((X<A)) == '-A__xy*X__y*e_x + A__xy*X__x*e_y'
        assert str((X>A)) == 'A*X__x*e_x + A*X__y*e_y'

        assert str((A*X)) == '(A*X__x + A__xy*X__y)*e_x + (A*X__y - A__xy*X__x)*e_y'
        assert str((A|X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'
        assert str((A<X)) == 'A*X__x*e_x + A*X__y*e_y'
        assert str((A>X)) == 'A__xy*X__y*e_x - A__xy*X__x*e_y'

        return

    def test_check_generalized_BAC_CAB_formulas(self):

        (a,b,c,d,e) = Ga('a b c d e').mv()

        assert str(a|(b*c)) == '-(a.c)*b + (a.b)*c'
        assert str(a|(b^c)) == '-(a.c)*b + (a.b)*c'
        assert str(a|(b^c^d)) == '(a.d)*b^c - (a.c)*b^d + (a.b)*c^d'

        expr = (a|(b^c))+(c|(a^b))+(b|(c^a)) # = (a.b)*c - (b.c)*a - ((a.b)*c - (b.c)*a)
        assert str(expr.simplify()) == ZERO_STR

        assert str(a*(b^c)-b*(a^c)+c*(a^b)) == '3*a^b^c'
        assert str(a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)) == '4*a^b^c^d'
        assert str((a^b)|(c^d)) == '-(a.c)*(b.d) + (a.d)*(b.c)'
        assert str(((a^b)|c)|d) == '-(a.c)*(b.d) + (a.d)*(b.c)'
        assert str(Ga.com(a^b,c^d)) == '-(b.d)*a^c + (b.c)*a^d + (a.d)*b^c - (a.c)*b^d'
        assert str((a|(b^c))|(d^e)) == '(-(a.b)*(c.e) + (a.c)*(b.e))*d + ((a.b)*(c.d) - (a.c)*(b.d))*e'

        return

    def test_derivatives_in_rectangular_coordinates(self):

        X = (x, y, z) = symbols('x y z')
        o3d = Ga('e_x e_y e_z', g=[1, 1, 1], coords=X)
        (ex, ey, ez) = o3d.mv()
        grad = o3d.grad

        f = o3d.mv('f', 'scalar', f=True)
        A = o3d.mv('A', 'vector', f=True)
        B = o3d.mv('B', 'bivector', f=True)
        C = o3d.mv('C', 'mv', f=True)

        assert str(f) == 'f'
        assert str(A) == 'A__x*e_x + A__y*e_y + A__z*e_z'
        assert str(B) == 'B__xy*e_x^e_y + B__xz*e_x^e_z + B__yz*e_y^e_z'
        assert str(C) == 'C + C__x*e_x + C__y*e_y + C__z*e_z + C__xy*e_x^e_y + C__xz*e_x^e_z + C__yz*e_y^e_z + C__xyz*e_x^e_y^e_z'

        assert str(grad*f) == 'D{x}f*e_x + D{y}f*e_y + D{z}f*e_z'
        assert str(grad|A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
        assert str(grad*A) == 'D{x}A__x + D{y}A__y + D{z}A__z + (-D{y}A__x + D{x}A__y)*e_x^e_y + (-D{z}A__x + D{x}A__z)*e_x^e_z + (-D{z}A__y + D{y}A__z)*e_y^e_z'

        assert str(-o3d.I()*(grad^A)) == '(-D{z}A__y + D{y}A__z)*e_x + (D{z}A__x - D{x}A__z)*e_y + (-D{y}A__x + D{x}A__y)*e_z'
        assert str(grad*B) == '(-D{y}B__xy - D{z}B__xz)*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z + (D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
        assert str(grad^B) == '(D{z}B__xy - D{y}B__xz + D{x}B__yz)*e_x^e_y^e_z'
        assert str(grad|B) == '(-D{y}B__xy - D{z}B__xz)*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'

        assert str(grad<A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
        assert str(grad>A) == 'D{x}A__x + D{y}A__y + D{z}A__z'
        assert str(grad<B) == '(-D{y}B__xy - D{z}B__xz)*e_x + (D{x}B__xy - D{z}B__yz)*e_y + (D{x}B__xz + D{y}B__yz)*e_z'
        assert str(grad>B) == ZERO_STR
        assert str(grad<C) == 'D{x}C__x + D{y}C__y + D{z}C__z + (-D{y}C__xy - D{z}C__xz)*e_x + (D{x}C__xy - D{z}C__yz)*e_y + (D{x}C__xz + D{y}C__yz)*e_z + D{z}C__xyz*e_x^e_y - D{y}C__xyz*e_x^e_z + D{x}C__xyz*e_y^e_z'
        assert str(grad>C) == 'D{x}C__x + D{y}C__y + D{z}C__z + D{x}C*e_x + D{y}C*e_y + D{z}C*e_z'

        return

    def test_derivatives_in_spherical_coordinates(self):

        X = (r, th, phi) = symbols('r theta phi')
        s3d = Ga('e_r e_theta e_phi', g=[1, r ** 2, r ** 2 * sin(th) ** 2], coords=X, norm=True)
        (er, eth, ephi) = s3d.mv()
        grad = s3d.grad

        f = s3d.mv('f', 'scalar', f=True)
        A = s3d.mv('A', 'vector', f=True)
        B = s3d.mv('B', 'bivector', f=True)

        assert str(f) == 'f'
        assert str(A) == 'A__r*e_r + A__theta*e_theta + A__phi*e_phi'
        assert str(B) == 'B__rtheta*e_r^e_theta + B__rphi*e_r^e_phi + B__thetaphi*e_theta^e_phi'

        assert str(grad*f) == 'D{r}f*e_r + D{theta}f*e_theta/r + D{phi}f*e_phi/(r*sin(theta))'
        assert str((grad|A).simplify()) == '(r*D{r}A__r + 2*A__r + A__theta/tan(theta) + D{theta}A__theta + D{phi}A__phi/sin(theta))/r'
        assert str(-s3d.I()*(grad^A)) == '(A__phi/tan(theta) + D{theta}A__phi - D{phi}A__theta/sin(theta))*e_r/r + (-r*D{r}A__phi - A__phi + D{phi}A__r/sin(theta))*e_theta/r + (r*D{r}A__theta + A__theta - D{theta}A__r)*e_phi/r'

        assert latex(grad) == r'\boldsymbol{e}_{r} \frac{\partial}{\partial r} + \boldsymbol{e}_{\theta } \frac{1}{r} \frac{\partial}{\partial \theta } + \boldsymbol{e}_{\phi } \frac{1}{r \sin{\left (\theta  \right )}} \frac{\partial}{\partial \phi }'
        assert latex(B|(eth^ephi)) == r'- B^{\theta \phi } {\left (r,\theta ,\phi  \right )}'

        assert str(grad^B) == '(r*D{r}B__thetaphi - B__rphi/tan(theta) + 2*B__thetaphi - D{theta}B__rphi + D{phi}B__rtheta/sin(theta))*e_r^e_theta^e_phi/r'

        return

    def test_rounding_numerical_components(self):

        o3d = Ga('e_x e_y e_z', g=[1, 1, 1])
        (ex, ey, ez) = o3d.mv()

        X = 1.2*ex+2.34*ey+0.555*ez
        Y = 0.333*ex+4*ey+5.3*ez

        assert str(X) == '1.2*e_x + 2.34*e_y + 0.555*e_z'
        assert str(Nga(X,2)) == '1.2*e_x + 2.3*e_y + 0.55*e_z'
        assert str(X*Y) == '12.7011000000000 + 4.02078*e_x^e_y + 6.175185*e_x^e_z + 10.182*e_y^e_z'
        assert str(Nga(X*Y,2)) == '13. + 4.0*e_x^e_y + 6.2*e_x^e_z + 10.0*e_y^e_z'

        return

    def test_noneuclidian_distance_calculation(self):
        from sympy import solve,sqrt

        g = '0 # #,# 0 #,# # 1'
        necl = Ga('X Y e',g=g)
        (X,Y,e) = necl.mv()

        assert str((X^Y)*(X^Y)) == '(X.Y)**2'

        L = X^Y^e
        B = (L*e).expand().blade_rep() # D&L 10.152
        assert str(B) == 'X^Y - (Y.e)*X^e + (X.e)*Y^e'
        Bsq = B*B
        assert str(Bsq) == '(X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e)'
        Bsq = Bsq.scalar()
        assert str(B) == 'X^Y - (Y.e)*X^e + (X.e)*Y^e'

        BeBr = B*e*B.rev()
        assert str(BeBr) == '(X.Y)*(-(X.Y) + 2*(X.e)*(Y.e))*e'
        assert str(B*B) == '(X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e)'
        assert str(L*L) == '(X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e)' # D&L 10.153

        (s,c,Binv,M,S,C,alpha) = symbols('s c (1/B) M S C alpha')

        XdotY = necl.g[0,1]
        Xdote = necl.g[0,2]
        Ydote = necl.g[1,2]

        Bhat = Binv*B # D&L 10.154
        R = c+s*Bhat # Rotor R = exp(alpha*Bhat/2)
        assert str(R) == 'c + (1/B)*s*X^Y - (1/B)*(Y.e)*s*X^e + (1/B)*(X.e)*s*Y^e'

        Z = R*X*R.rev() # D&L 10.155
        Z.obj = expand(Z.obj)
        Z.obj = Z.obj.collect([Binv,s,c,XdotY])
        assert str(Z) == '((1/B)**2*(X.Y)**2*s**2 - 2*(1/B)**2*(X.Y)*(X.e)*(Y.e)*s**2 + 2*(1/B)*(X.Y)*c*s - 2*(1/B)*(X.e)*(Y.e)*c*s + c**2)*X + 2*(1/B)*(X.e)**2*c*s*Y + 2*(1/B)*(X.Y)*(X.e)*s*(-(1/B)*(X.Y)*s + 2*(1/B)*(X.e)*(Y.e)*s - c)*e'
        W = Z|Y
        # From this point forward all calculations are with sympy scalars
        W = W.scalar()
        assert str(W) == '(1/B)**2*(X.Y)**3*s**2 - 4*(1/B)**2*(X.Y)**2*(X.e)*(Y.e)*s**2 + 4*(1/B)**2*(X.Y)*(X.e)**2*(Y.e)**2*s**2 + 2*(1/B)*(X.Y)**2*c*s - 4*(1/B)*(X.Y)*(X.e)*(Y.e)*c*s + (X.Y)*c**2'
        W = expand(W)
        W = simplify(W)
        W = W.collect([s*Binv])

        M = 1/Bsq
        W = W.subs(Binv**2,M)
        W = simplify(W)
        Bmag = sqrt(XdotY**2-2*XdotY*Xdote*Ydote)
        W = W.collect([Binv*c*s,XdotY])

        #Double angle substitutions

        W = W.subs(2*XdotY**2-4*XdotY*Xdote*Ydote,2/(Binv**2))
        W = W.subs(2*c*s,S)
        W = W.subs(c**2,(C+1)/2)
        W = W.subs(s**2,(C-1)/2)
        W = simplify(W)
        W = W.subs(Binv,1/Bmag)
        W = expand(W)

        assert str(W.simplify()) == '(X.Y)*C - (X.e)*(Y.e)*C + (X.e)*(Y.e) + S*sqrt((X.Y)*((X.Y) - 2*(X.e)*(Y.e)))'

        Wd = collect(W,[C,S],exact=True,evaluate=False)

        Wd_1 = Wd[one]
        Wd_C = Wd[C]
        Wd_S = Wd[S]

        assert str(Wd_1) == '(X.e)*(Y.e)'
        assert str(Wd_C) == '(X.Y) - (X.e)*(Y.e)'

        assert str(Wd_S) == 'sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'
        assert str(Bmag) == 'sqrt((X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e))'

        Wd_1 = Wd_1.subs(Binv,1/Bmag)
        Wd_C = Wd_C.subs(Binv,1/Bmag)
        Wd_S = Wd_S.subs(Binv,1/Bmag)

        lhs = Wd_1+Wd_C*C
        rhs = -Wd_S*S
        lhs = lhs**2
        rhs = rhs**2
        W = expand(lhs-rhs)
        W = expand(W.subs(1/Binv**2,Bmag**2))
        W = expand(W.subs(S**2,C**2-1))
        W = W.collect([C,C**2],evaluate=False)

        a = simplify(W[C**2])
        b = simplify(W[C])
        c = simplify(W[one])

        assert str(a) == '(X.e)**2*(Y.e)**2'
        assert str(b) == '2*(X.e)*(Y.e)*((X.Y) - (X.e)*(Y.e))'
        assert str(c) == '(X.Y)**2 - 2*(X.Y)*(X.e)*(Y.e) + (X.e)**2*(Y.e)**2'

        x = Symbol('x')
        C =  solve(a*x**2+b*x+c,x)[0]
        assert str(expand(simplify(expand(C)))) == '-(X.Y)/((X.e)*(Y.e)) + 1'

        return

    def test_conformal_representations_of_circles_lines_spheres_and_planes(self):
        global n,nbar

        g = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

        cnfml3d = Ga('e_1 e_2 e_3 n nbar',g=g)

        (e1,e2,e3,n,nbar) = cnfml3d.mv()

        e = n+nbar

        #conformal representation of points

        A = make_vector(e1,ga=cnfml3d)    # point a = (1,0,0)  A = F(a)
        B = make_vector(e2,ga=cnfml3d)    # point b = (0,1,0)  B = F(b)
        C = make_vector(-e1,ga=cnfml3d)   # point c = (-1,0,0) C = F(c)
        D = make_vector(e3,ga=cnfml3d)    # point d = (0,0,1)  D = F(d)
        X = make_vector('x',3,ga=cnfml3d)

        assert str(A) == 'e_1 + n/2 - nbar/2'
        assert str(B) == 'e_2 + n/2 - nbar/2'
        assert str(C) == '-e_1 + n/2 - nbar/2'
        assert str(D) == 'e_3 + n/2 - nbar/2'
        assert str(X) == 'x1*e_1 + x2*e_2 + x3*e_3 + (x1**2/2 + x2**2/2 + x3**2/2)*n - nbar/2'

        assert str((A^B^C^X)) == '-x3*e_1^e_2^e_3^n + x3*e_1^e_2^e_3^nbar + (x1**2/2 + x2**2/2 + x3**2/2 - 1/2)*e_1^e_2^n^nbar'
        assert str((A^B^n^X)) == '-x3*e_1^e_2^e_3^n + (x1/2 + x2/2 - 1/2)*e_1^e_2^n^nbar + x3*e_1^e_3^n^nbar/2 - x3*e_2^e_3^n^nbar/2'
        assert str((((A^B)^C)^D)^X) == '(-x1**2/2 - x2**2/2 - x3**2/2 + 1/2)*e_1^e_2^e_3^n^nbar'
        assert str((A^B^n^D^X)) == '(-x1/2 - x2/2 - x3/2 + 1/2)*e_1^e_2^e_3^n^nbar'

        L = (A^B^e)^X

        assert str(L) == '-x3*e_1^e_2^e_3^n - x3*e_1^e_2^e_3^nbar + (-x1**2/2 + x1 - x2**2/2 + x2 - x3**2/2 - 1/2)*e_1^e_2^n^nbar + x3*e_1^e_3^n^nbar - x3*e_2^e_3^n^nbar'

        return

    def test_properties_of_geometric_objects(self):

        global n, nbar

        g = '# # # 0 0,'+ \
            '# # # 0 0,'+ \
            '# # # 0 0,'+ \
            '0 0 0 0 2,'+ \
            '0 0 0 2 0'

        c3d = Ga('p1 p2 p3 n nbar',g=g)

        (p1,p2,p3,n,nbar) = c3d.mv()

        P1 = F(p1)
        P2 = F(p2)
        P3 = F(p3)

        L = P1^P2^n
        delta = (L|n)|nbar
        assert str(delta) == '2*p1 - 2*p2'

        C = P1^P2^P3
        delta = ((C^n)|n)|nbar
        assert str(delta) == '2*p1^p2 - 2*p1^p3 + 2*p2^p3'
        assert str((p2-p1)^(p3-p1)) == 'p1^p2 - p1^p3 + p2^p3'

        return

    def test_extracting_vectors_from_conformal_2_blade(self):

        g = '0 -1 #,'+ \
            '-1 0 #,'+ \
            '# # #'

        e2b = Ga('P1 P2 a',g=g)

        (P1,P2,a) = e2b.mv()

        B = P1^P2
        Bsq = B*B
        assert str(Bsq) == '1'
        ap = a-(a^B)*B
        assert str(ap) == '-(P2.a)*P1 - (P1.a)*P2'

        Ap = ap+ap*B
        Am = ap-ap*B

        assert str(Ap) == '-2*(P2.a)*P1'
        assert str(Am) == '-2*(P1.a)*P2'

        assert str(Ap*Ap) == ZERO_STR
        assert str(Am*Am) == ZERO_STR

        aB = a|B
        assert str(aB) == '-(P2.a)*P1 + (P1.a)*P2'

        return

    def test_reciprocal_frame_test(self):

        g = '1 # #,'+ \
            '# 1 #,'+ \
            '# # 1'

        g3dn = Ga('e1 e2 e3',g=g)

        (e1,e2,e3) = g3dn.mv()

        E = e1^e2^e3
        Esq = (E*E).scalar()
        assert str(E) == 'e1^e2^e3'
        assert str(Esq) == '(e1.e2)**2 - 2*(e1.e2)*(e1.e3)*(e2.e3) + (e1.e3)**2 + (e2.e3)**2 - 1'
        Esq_inv = 1/Esq

        E1 = (e2^e3)*E
        E2 = (-1)*(e1^e3)*E
        E3 = (e1^e2)*E

        assert str(E1) == '((e2.e3)**2 - 1)*e1 + ((e1.e2) - (e1.e3)*(e2.e3))*e2 + (-(e1.e2)*(e2.e3) + (e1.e3))*e3'
        assert str(E2) == '((e1.e2) - (e1.e3)*(e2.e3))*e1 + ((e1.e3)**2 - 1)*e2 + (-(e1.e2)*(e1.e3) + (e2.e3))*e3'
        assert str(E3) == '(-(e1.e2)*(e2.e3) + (e1.e3))*e1 + (-(e1.e2)*(e1.e3) + (e2.e3))*e2 + ((e1.e2)**2 - 1)*e3'

        w = (E1|e2)
        w = w.expand()
        assert str(w) == ZERO_STR

        w = (E1|e3)
        w = w.expand()
        assert str(w) == ZERO_STR

        w = (E2|e1)
        w = w.expand()
        assert str(w) == ZERO_STR

        w = (E2|e3)
        w = w.expand()
        assert str(w) == ZERO_STR

        w = (E3|e1)
        w = w.expand()
        assert str(w) == ZERO_STR

        w = (E3|e2)
        w = w.expand()
        assert str(w) == ZERO_STR

        w = (E1|e1)
        w = (w.expand()).scalar()
        Esq = expand(Esq)
        assert str(simplify(w/Esq)) == '1'

        w = (E2|e2)
        w = (w.expand()).scalar()
        assert str(simplify(w/Esq)) == '1'

        w = (E3|e3)
        w = (w.expand()).scalar()
        assert str(simplify(w/Esq)) == '1'

        return
