from __future__ import print_function
from galgebra.printer import xpdf,Get_Program,Print_Function,Format
from galgebra.ga import Ga

def basic_multivector_operations_3D():
    Print_Function()

    (g3d,ex,ey,ez) = Ga.build('e*x|y|z')

    print('g_{ij} =',g3d.g)

    A = g3d.mv('A','mv')

    print(A.Fmt(1,'A'))
    print(A.Fmt(2,'A'))
    print(A.Fmt(3,'A'))

    print(A.even().Fmt(1,'%A_{+}'))
    print(A.odd().Fmt(1,'%A_{-}'))

    X = g3d.mv('X','vector')
    Y = g3d.mv('Y','vector')

    print(X.Fmt(1,'X'))
    print(Y.Fmt(1,'Y'))

    print((X*Y).Fmt(2,'X*Y'))
    print((X^Y).Fmt(2,'X^Y'))
    print((X|Y).Fmt(2,'X|Y'))
    return

def basic_multivector_operations_2D():
    Print_Function()

    (g2d,ex,ey) = Ga.build('e*x|y')

    print('g_{ij} =',g2d.g)

    X = g2d.mv('X','vector')
    A = g2d.mv('A','spinor')

    print(X.Fmt(1,'X'))
    print(A.Fmt(1,'A'))

    print((X|A).Fmt(2,'X|A'))
    print((X<A).Fmt(2,'X<A'))
    print((A>X).Fmt(2,'A>X'))
    return

def dummy():
    return

def main():
    Get_Program(True)
    Format()

    basic_multivector_operations_3D()
    basic_multivector_operations_2D()

    # xpdf('simple_test_latex.tex')
    xpdf('simple_check_latex.tex', pdfprog=None)
    return

if __name__ == "__main__":
    main()
