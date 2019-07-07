from __future__ import print_function
from sympy import *
from galgebra.printer import Format,xpdf,xdvi
from galgebra.ga import Ga

def main():
    Format()
    snr=1
    g = '0 0 1 0 ,0 0 0 1 ,1 0 0 0 ,0 1 0 0'
    sk4coords = (e1,e2,e3,e4) = symbols('e1 e2 e3 e4')
    sk4 = Ga('e_1 e_2 e_3 e_4', g=g, coords=sk4coords)
    (e1,e2,e3,e4) = sk4.mv()
    print('g_{ii} =',sk4.g)

    v = symbols('v', real=True)
    x1=(e1+e3)/sqrt(2)
    x2=(e2+e4)/sqrt(2)
    print('x_1<x_1==',x1<x1)
    print('x_1<x_2==',x1<x2)
    print('x_2<x_1==',x2<x1)
    print('x_2<x_2==',x2<x2)
    print(r'#$-\infty < v < \infty$')
    print('(-v*(x_1^x_2)/2).exp()==',(-v*(x1^x2)/2).exp())
    v = symbols('v', real=True, positive=True)
    print(r'#$0\le v < \infty$')
    print('(-v*(x_1^x_2)/2).exp()==',(-v*(x1^x2)/2).exp())

    xpdf(pdfprog=None)
    return

if __name__ == "__main__":
    main()
