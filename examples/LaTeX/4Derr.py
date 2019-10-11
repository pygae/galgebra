from __future__ import print_function
from sympy import *
from galgebra.printer import Format,xtex
from galgebra.ga import Ga

def main():
    Format()
    snr=1
    g = '0 0 1 0 ,0 0 0 1 ,1 0 0 0 ,0 1 0 0'
    sk4coords = (e1,e2,e3,e4) = symbols('e1 e2 e3 e4', real=True)
    sk4 = Ga('e_1 e_2 e_3 e_4', g=g, coords=sk4coords)
    (e1,e2,e3,e4) = sk4.mv()
    print('g =',sk4.g)

    v = symbols('v', real=True)
    x1=(e1+e3)/sqrt(2)
    x2=(e2+e4)/sqrt(2)
    print(r'x_1\lfloor x_1==',x1<x1)
    print('x_1\lfloor x_2==',x1<x2)
    print('x_2\lfloor x_1==',x2<x1)
    print('x_2\lfloor x_2==',x2<x2)
    print(r'-\infty < v < \infty')
    print(r'e^{-v(x_1\W x_2)/2} ==',(-v*(x1^x2)/2).exp())
    v = symbols('v', real=True, positive=True)
    print(r'0\le v \le \infty')
    print(r'e^{-v(x_1\W x_2)/2)} ==',(-v*(x1^x2)/2).exp())

    xtex()
    return

if __name__ == "__main__":
    main()
