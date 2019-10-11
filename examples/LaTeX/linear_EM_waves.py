from __future__ import print_function
import sys
from sympy import symbols,sin,cos,exp,I,Matrix,solve,simplify
from galgebra.printer import Format,xtex,Get_Program,Print_Function,hline
from galgebra.ga import Ga
from galgebra.metric import linear_expand

def EM_Waves_in_Geom_Calculus_Complex():
    #Print_Function()
    X = (t,x,y,z) = symbols('t x y z',real=True)
    g = '1 # # 0,# 1 # 0,# # 1 0,0 0 0 -1'
    coords = (xE,xB,xk,t) = symbols('x_E x_B x_k t',real=True)
    (EBkst,eE,eB,ek,et) = Ga.build('e_E e_B e_k e_t',g=g,coords=coords)

    i = EBkst.i

    E,B,k,w = symbols('E B k omega',real=True)

    F = E*eE*et+i*B*eB*et
    K = k*ek+w*et
    X = xE*eE+xB*eB+xk*ek+t*et
    KX = (K|X).scalar()
    F = F*exp(I*KX)

    g = EBkst.g

    print('g =', g)
    print('X =', X)
    print('K =', K)
    print('K|X =', KX)
    hline()
    print('F =', F.Fmt(3))
    hline()
    gradF = EBkst.grad*F

    gradF = gradF.simplify()

    print(r'\nabla F =', gradF.Fmt(3))
    hline()
    gradF = gradF.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    KX = KX.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    print(r'\mbox{Substituting }e_{E}\cdot e_{B} = e_{E}\cdot e_{k} = e_{B}\cdot e_{k} = 0')

    print(r'\lp\bm{\nabla}F\rp/\lp ie^{iK\cdot X}\rp = ',(gradF / (I*exp(I*KX))).Fmt(3))
    hline()
    return

def EM_Waves_in_Geom_Calculus_Real():
    #Print_Function()
    X = (t,x,y,z) = symbols('t x y z',real=True)
    g = '1 # # 0,# 1 # 0,# # 1 0,0 0 0 -1'
    coords = (xE,xB,xk,t) = symbols('x_E x_B x_k t',real=True)
    (EBkst,eE,eB,ek,et) = Ga.build('e_E e_B e_k e_t',g=g,coords=coords)

    i = EBkst.i

    E,B,k,w = symbols('E B k omega',real=True)

    F = E*eE*et+i*B*eB*et
    K = k*ek+w*et
    X = xE*eE+xB*eB+xk*ek+t*et
    KX = (K|X).scalar()
    F = F*sin(KX)

    g = EBkst.g

    print('g =', g)
    print('X =', X)
    print('K =', K)
    print('K|X =', KX)
    hline()
    print('F =',F.Fmt(3))
    hline()
    gradF = EBkst.grad*F

    gradF = gradF.simplify()

    print(r'\nabla F =',(gradF).Fmt(3))
    hline()
    gradF = gradF.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    KX = KX.subs({g[0,1]:0,g[0,2]:0,g[1,2]:0})

    print(r'\mbox{Substituting }e_{E}\cdot e_{B} = e_{E}\cdot e_{k} = e_{B}\cdot e_{k} = 0')

    print(r'\lp\bm{\nabla}F\rp/\lp \f{\cos}{K\cdot X}\rp =',(gradF / (cos(KX))).Fmt(3))
    hline()
    return

def dummy():
    return

def main():
    Get_Program()
    Format()

    EM_Waves_in_Geom_Calculus_Complex()
    EM_Waves_in_Geom_Calculus_Real()
    xtex()
    return

if __name__ == "__main__":
    main()
