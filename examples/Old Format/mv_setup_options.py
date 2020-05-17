from __future__ import print_function
from sympy import symbols
from galgebra.deprecated import MV
from galgebra.printer import enhance_print,Print_Function

def MV_setup_options():
    (e1,e2,e3) = MV.setup('e_1 e_2 e_3','[1,1,1]')
    v = MV('v', 'vector')
    print(v)

    (e1,e2,e3) = MV.setup('e*1|2|3','[1,1,1]')
    v = MV('v', 'vector')
    print(v)

    (e1,e2,e3) = MV.setup('e*x|y|z','[1,1,1]')
    v = MV('v', 'vector')
    print(v)

    coords = symbols('x y z')
    (e1,e2,e3,grad) = MV.setup('e','[1,1,1]',coords=coords)
    v = MV('v', 'vector')
    print(v)

    return

def main():
    enhance_print()
    MV_setup_options()
    return

if __name__ == "__main__":
    main()
