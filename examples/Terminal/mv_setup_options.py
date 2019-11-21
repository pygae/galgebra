from __future__ import print_function
from sympy import symbols
from galgebra.ga import Ga
from galgebra.printer import Eprint, Get_Program, Print_Function

def Mv_setup_options():
    Print_Function()

    (o3d,e1,e2,e3) = Ga.build('e_1 e_2 e_3',g=[1,1,1])
    v = o3d.mv('v', 'vector')
    print(v)

    (o3d,e1,e2,e3) = Ga.build('e*1|2|3',g=[1,1,1])
    v = o3d.mv('v', 'vector')
    print(v)

    (o3d,e1,e2,e3) = Ga.build('e*x|y|z',g=[1,1,1])
    v = o3d.mv('v', 'vector')
    print(v)

    coords = symbols('x y z',real=True)
    (o3d,e1,e2,e3) = Ga.build('e',g=[1,1,1],coords=coords)
    v = o3d.mv('v', 'vector')
    print(v)

    print(v.grade(2))
    print(v.i_grade)

    return

def dummy():
    return

def main():
    Get_Program()
    Eprint()
    Mv_setup_options()
    return

if __name__ == "__main__":
    main()
