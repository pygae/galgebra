from __future__ import print_function

import sys
from sympy import symbols,sin,cos,exp,I,Matrix,simplify,Eq,S
from sympy.solvers import solve
from galgebra.printer import Format,xtex,Get_Program,Print_Function,Fmt
from galgebra.ga import Ga
from galgebra.metric import linear_expand

def EM_Waves_in_Geom_Calculus():
    #Print_Function()
    X = (t,x,y,z) = symbols('t x y z',real=True)
    (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=X)

    i = st4d.I()

    B = st4d.mv('B','vector')
    E = st4d.mv('E','vector')
    B.set_coef(1,0,0)
    E.set_coef(1,0,0)
    B *= g0
    E *= g0
    F = E+i*B

    kx, ky, kz, w = symbols('k_x k_y k_z omega',real=True)
    kv = kx*g1+ky*g2+kz*g3
    xv = x*g1+y*g2+z*g3
    KX = ((w*g0+kv)|(t*g0+xv)).scalar()

    Ixyz = g1*g2*g3

    F = F*exp(I*KX)

    print(r'\T{Pseudo Scalar\;\;}I =',i)
    print(r'I_{xyz} =',Ixyz)
    print(r'\T{Electromagnetic Field Bi-Vector\;\;} F =',F.Fmt(3))
    gradF = st4d.grad*F

    print(r'\T{Geom Derivative of Electomagnetic Field Bi-Vector}')
    print(r'\nabla F = 0 =',gradF.Fmt(3))

    gradF = gradF / (I * exp(I*KX))
    print(r'\lp\bs{\nabla}F\rp /\lp i e^{iK\cdot X}\rp = 0 =',gradF.Fmt(3))

    g = '-1 # 0 0,# -1 0 0,0 0 -1 0,0 0 0 1'
    X = (xE,xB,xk,t) = symbols('x_E x_B x_k t',real=True)

    # The correct signature sig=p for signature (p,q) is needed
    # since the correct value of I**2 is needed to compute exp(I)
    (EBkst,eE,eB,ek,et) = Ga.build('e_E e_B e_k t',g=g,coords=X,sig=1)

    i = EBkst.I()

    E,B,k,w = symbols('E B k omega',real=True)

    F = E*eE*et+i*B*eB*et
    kv = k*ek+w*et
    xv = xE*eE+xB*eB+xk*ek+t*et
    KX = (kv|xv).scalar()
    F = F*exp(I*KX)

    print(r'\mbox{set } e_{E}\cdot e_{k} = e_{B}\cdot e_{k} = 0'+\
           r'\mbox{ and } e_{E}\cdot e_{E} = e_{B}\cdot e_{B} = e_{k}\cdot e_{k} = -e_{t}\cdot e_{t} = 1')

    print('g =', EBkst.g)

    print('K|X =',KX)
    print('F =',F)
    print(r'\nabla F = 0 =',(EBkst.grad*F).Fmt(3))

    gradF_reduced = (EBkst.grad*F)/(I*exp(I*KX))

    print(r'\lp\bs{\nabla}F\rp/\lp ie^{iK\cdot X} \rp = 0 =',gradF_reduced.Fmt(3))

    print(r'\mbox{Previous equation requires that: }e_{E}\cdot e_{B} = 0\mbox{ if }B\ne 0\mbox{ and }k\ne 0')

    gradF_reduced = gradF_reduced.subs({EBkst.g[0,1]:0}).expand()
    print(r'\lp\bs{\nabla}F\rp/\lp ie^{iK\cdot X} \rp = 0 =',gradF_reduced.Fmt(3))

    (coefs,bases) = linear_expand(gradF_reduced.obj)

    eq1 = coefs[0]
    eq2 = coefs[1]

    print(r'0 =', eq1)
    print(r'0 =', eq2)

    B1 = solve(eq1,B)[0]
    B2 = solve(eq2,B)[0]

    eq3 = B1-B2

    print(r'\mbox{eq3 = eq1-eq2: }0 =',eq3)
    eq3 = simplify(eq3 / E)
    print(r'\mbox{eq3 = (eq1-eq2)/E: }0 =',eq3)
    print('k =',Matrix(solve(eq3,k)))
    print('B =',Matrix([B1.subs(w,k),B1.subs(-w,k)]))

    return

def dummy():
    return

def main():
    #Get_Program()
    Format()

    EM_Waves_in_Geom_Calculus()

    # xpdf()
    xtex()
    return

if __name__ == "__main__":
    main()
