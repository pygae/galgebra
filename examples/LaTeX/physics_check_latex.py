from __future__ import print_function

import sys
from sympy import symbols,sin,cos
from galgebra.printer import Format,xtex,Get_Program,Print_Function
from galgebra.ga import Ga

def Maxwells_Equations_in_Geom_Calculus():
    Print_Function()
    X = symbols('t x y z',real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=X)

    I = st4d.i

    B = st4d.mv('B','vector',f=True)
    E = st4d.mv('E','vector',f=True)
    B.set_coef(1,0,0)
    E.set_coef(1,0,0)
    B *= g0
    E *= g0
    J = st4d.mv('J','vector',f=True)
    F = E+I*B

    print(r'\T{Pseudo Scalar\;\;}I =',I)
    print(r'\T{Magnetic Field Bi-Vector\;\;} B = \bs{B\gamma_{t}} =',B)
    print(r'\T{Electric Field Bi-Vector\;\;} E = \bs{E\gamma_{t}} =',E)
    print(r'\T{Electromagnetic Field Bi-Vector\;\;} F = E+IB =',F)
    print(r'\T{Four Current Density\;\;} J =',J)
    gradF = st4d.grad*F
    print(r'\T{Geom Derivative of Electomagnetic Field Bi-Vector}')
    print(r'\nabla F =',gradF.Fmt(3))

    print(r'\T{Maxwell Equations}')
    print(r'\nabla F = J')
    print(r'\T{Div $E$ and Curl $H$ Equations}')
    print(r'\grade{\nabla F}{1} -J = 0 =',(gradF.get_grade(1)-J).Fmt(3))
    print(r'\T{Curl $E$ and Div $B$ equations}')
    print(r'\grade{\nabla F}{3} = 0 =',(gradF.get_grade(3)).Fmt(3))
    return

def Dirac_Equation_in_Geom_Calculus():
    Print_Function()
    coords = symbols('t x y z',real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)
    I = st4d.i

    (m,e) = symbols('m e')

    psi = st4d.mv('psi','spinor',f=True)
    A = st4d.mv('A','vector',f=True)
    sig_z = g3*g0

    print(r'\T{4-Vector Potential\;\;}\\bs{A} =',A)
    print(r'\T{8-component real spinor\;\;}\bs{\psi} =',psi)

    dirac_eq = (st4d.grad*psi)*I*sig_z-e*A*psi-m*psi*g0
    dirac_eq = dirac_eq.simplify()

    print(r'\T{Dirac Equation\;\;}\nabla \bs{\psi} I \sigma_{z}-e\bs{A}\bs{\psi}-m\bs{\psi}\gamma_{t} = 0 =',dirac_eq.Fmt(3))

    return

def Lorentz_Tranformation_in_Geog_Algebra():
    Print_Function()
    (alpha,beta,gamma) = symbols('alpha beta gamma')
    (x,t,xp,tp) = symbols("x t x' t'",real=True)
    (st2d,g0,g1) = Ga.build('gamma*t|x',g=[1,-1])

    from sympy import sinh,cosh

    R = cosh(alpha/2)+sinh(alpha/2)*(g0^g1)
    X = t*g0+x*g1
    Xp = tp*g0+xp*g1
    print('R =',R)

    print(r"t\bm{\gamma_{t}}+x\bm{\gamma_{x}} = t'\bm{\gamma'_{t}}+x'\bm{\gamma'_{x}} = R\lp t'\bm{\gamma_{t}}+x'\bm{\gamma_{x}}\rp R^{\dagger}")

    Xpp = R*Xp*R.rev()
    Xpp = Xpp.collect()
    Xpp = Xpp.trigsimp()
    print(r't\bm{\gamma_{t}}+x\bm{\gamma_{x}} =',Xpp)
    Xpp = Xpp.subs({sinh(alpha):gamma*beta,cosh(alpha):gamma})

    print(r'\f{\sinh}{\alpha} = \gamma\beta')
    print(r'\f{\cosh}{\alpha} = \gamma')

    print(r't\bm{\gamma_{t}}+x\bm{\gamma_{x}} =',Xpp.collect())
    return

def General_Lorentz_Tranformation():
    Print_Function()
    (alpha,beta,gamma) = symbols('alpha beta gamma')
    (x,y,z,t) = symbols("x y z t",real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1])

    B = (x*g1+y*g2+z*g3)^(t*g0)
    print(B)
    print(B.exp(hint='+'))
    print(B.exp(hint='-'))

def Lie_Group():
    Print_Function()
    coords = symbols('t x y z',real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)
    I = st4d.i

    a = st4d.mv('a','vector')
    B = st4d.mv('B','bivector')
    print('a =',a)
    print('B =',B)
    print(r'a\cdot B =', a|B)
    print(r'(a\cdot B)\cdot B =',((a|B)|B).simplify().Fmt(3))
    print(r'((a\cdot B)\cdot B)\cdot B =',(((a|B)|B)|B).simplify().Fmt(3))

    return

def main():
    Get_Program()
    Format()

    Maxwells_Equations_in_Geom_Calculus()
    Dirac_Equation_in_Geom_Calculus()
    Lorentz_Tranformation_in_Geog_Algebra()
    General_Lorentz_Tranformation()
    Lie_Group()

    xtex()
    return

if __name__ == "__main__":
    main()
