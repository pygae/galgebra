from __future__ import print_function
import sys

from sympy import symbols,sin,cos,simplify
from galgebra.deprecated import MV
from galgebra.printer import enhance_print

def main():
    enhance_print()

    (ex,ey,ez) = MV.setup('e*x|y|z',metric='[1,1,1]')

    u = MV('u','vector')
    v = MV('v','vector')
    w = MV('w','vector')
    print(u)
    print(v)
    print(w)

    uv = u^v
    print(uv)
    print(uv.is_blade())
    uvw = u^v^w
    print(uvw)
    print(uvw.is_blade())

    print(simplify((uv*uv).scalar()))
    return

if __name__ == "__main__":
    main()
