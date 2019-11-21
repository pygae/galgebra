from __future__ import print_function
from sympy import symbols, Matrix
from galgebra.printer import xpdf, Format

def main():
    Format()
    a = Matrix ( 2, 2, ( 1, 2, 3, 4 ) )
    b = Matrix ( 2, 1, ( 5, 6 ) )
    c = a * b
    print(a,b,'=',c)

    x, y = symbols ( 'x, y' )

    d = Matrix ( 1, 2, ( x ** 3, y ** 3 ))
    e = Matrix ( 2, 2, ( x ** 2, 2 * x * y, 2 * x * y, y ** 2 ) )
    f = d * e

    print('%',d,e,'=',f)

    # xpdf()
    xpdf(pdfprog=None)
    return

if __name__ == "__main__":
    main()
