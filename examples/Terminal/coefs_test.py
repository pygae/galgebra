from __future__ import print_function
from sympy import symbols,S
from galgebra.ga import Ga
from galgebra.printer import Eprint, Get_Program, Print_Function

def coefs_test():
    Print_Function()

    (o3d,e1,e2,e3) = Ga.build('e_1 e_2 e_3',g=[1,1,1])
    print(o3d.blades.flat)
    print(o3d.mv_blades.flat)
    v = o3d.mv('v', 'vector')
    print(v)
    print(v.blade_coefs([e3,e1]))
    A = o3d.mv('A', 'mv')
    print(A)
    print(A.blade_coefs([e1^e3,e3,e1^e2,e1^e2^e3]))
    print(A.blade_coefs())
    return

def dummy():
    return

def main():
    Get_Program()
    Eprint()
    coefs_test()
    return

if __name__ == "__main__":
    main()
