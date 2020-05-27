from __future__ import print_function
from sympy import symbols
from galgebra.ga import Ga
from galgebra.printer import Format,xtex

def main():
    Format()

    coords = (x,y,z) = symbols('x y z',real=True)

    (o3d,ex,ey,ez) = Ga.build('e*x|y|z',g=[1,1,1],coords=coords)

    s = o3d.mv('s','scalar')
    v = o3d.mv('v','vector')
    b = o3d.mv('b','bivector')

    print(r'\T{3D Orthogonal Metric\newline}')

    print(r'\T{Multvectors:}')
    print('s =',s)
    print('v =',v)
    print('b =',b)

    print(r'\T{Products:}')

    X = ((s,'s'),(v,'v'),(b,'b'))

    for xi in X:
        print('')
        for yi in X:
            print(xi[1]+' '+yi[1]+' =',xi[0]*yi[0])
            print(xi[1]+r'\W '+yi[1]+' =',xi[0]^yi[0])
            if xi[1] != 's' and yi[1] != 's':
                print(xi[1]+r'\cdot '+yi[1]+' =',xi[0]|yi[0])
            print(xi[1]+r'\lfloor '+yi[1]+' =',xi[0]<yi[0])
            print(xi[1]+r'\rfloor '+yi[1]+' =',xi[0]>yi[0])

    fs = o3d.mv('s','scalar',f=True)
    fv = o3d.mv('v','vector',f=True)
    fb = o3d.mv('b','bivector',f=True)

    print(r'\T{Multivector Functions:}')

    print('s(X) =',fs)
    print('v(X) =',fv)
    print('b(X) =',fb)

    print(r'\T{Products:}')

    fX = ((o3d.grad,r'\nabla'),(fs,'s'),(fv,'v'),(fb,'b'))

    for xi in fX:
        print('')
        for yi in fX:
            if xi[1] == r'\nabla' and yi[1] == r'\nabla':
                pass
            else:
                print(xi[1]+' '+yi[1]+' =',xi[0]*yi[0])
                print(xi[1]+r'\W '+yi[1]+' =',xi[0]^yi[0])
                if xi[1] != 's' and yi[1] != 's':
                    print(xi[1]+r'\cdot '+yi[1]+' =',xi[0]|yi[0])
                print(xi[1]+r'\lfloor '+yi[1]+' =' ,xi[0]<yi[0])
                print(xi[1]+r'\rfloor '+yi[1]+' =' ,xi[0]>yi[0])


    (g2d,ex,ey) = Ga.build('e',coords=(x,y))

    print(r'\T{General 2D Metric\newline}')
    print(r'\T{Multivector Functions:}')

    s = g2d.mv('s','scalar',f=True)
    v = g2d.mv('v','vector',f=True)
    b = g2d.mv('v','bivector',f=True)

    print('s(X) =',s)
    print('v(X) =',v)
    print('b(X) =',b)

    X = ((g2d.grad,r'\nabla'),(s,'s'),(v,'v'))

    print(r'\T{Products:}')

    for xi in X:
        print('')
        for yi in X:
            if xi[1] == r'\nabla' and yi[1] == r'\nabla':
                pass
            else:
                print(xi[1]+' '+yi[1]+' =',xi[0]*yi[0])
                print(xi[1]+r'\W '+yi[1]+' =',xi[0]^yi[0])
                if xi[1] != 's' and yi[1] != 's':
                    print(xi[1]+r'\cdot '+yi[1]+' =',xi[0]|yi[0])
                else:
                    print(xi[1]+r'\cdot '+yi[1]+' = Not Allowed')
                print(xi[1]+r'\lfloor '+yi[1]+' =',xi[0]<yi[0])
                print(xi[1]+r'\rfloor '+yi[1]+' ='  ,xi[0]>yi[0])

    xtex()
    return

if __name__ == "__main__":
    main()
