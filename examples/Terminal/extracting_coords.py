"""
Extracting coordinates from multivectors
=========================================

Demonstrates the various ways to get scalar coefficients out of a
multivector, answering the question in issue 483.

Methods shown:
  - blade_coefs()  : coefficients for a list of basis blades
  - components()   : list of single-blade Mv objects
  - get_coefs(r)   : coefficients of all grade-r blades
  - scalar()       : the grade-0 part as a SymPy expression
  - get_grade(r)   : the grade-r part as an Mv

For numeric work (e.g. passing to math.sin), call float() or
sympy.N() on the SymPy expressions returned by blade_coefs().
"""

from sympy import symbols, cos, sin, pi
from galgebra.ga import Ga
from galgebra.printer import Eprint


def main():
    Eprint()

    # -- Setup: 3D Euclidean algebra --
    ga, e1, e2, e3 = Ga.build('e_1 e_2 e_3', g=[1, 1, 1])

    # A concrete vector
    v = 3 * e1 - 2 * e2 + 5 * e3
    print('v =', v)

    # 1) blade_coefs: get coefficients for specific blades
    print('\n--- blade_coefs ---')
    coefs = v.blade_coefs([e1, e2, e3])
    print('coefs of v along e1, e2, e3:', coefs)
    # For numeric use:
    print('float(coef of e1):', float(coefs[0]))

    # You can also request all blades (returns zeros for missing grades)
    A = ga.mv('A', 'mv')
    print('\nA =', A)
    print('all blade_coefs of A:', A.blade_coefs())

    # 2) components: split into single-blade multivectors
    print('\n--- components ---')
    for comp in v.components():
        print(' ', comp)

    # 3) get_coefs: coefficients of a specific grade
    print('\n--- get_coefs ---')
    B = 2 * (e1 ^ e2) - (e2 ^ e3)
    print('B =', B)
    print('grade-2 coefs of B:', B.get_coefs(2))

    # 4) scalar and get_grade
    print('\n--- scalar / get_grade ---')
    M = ga.mv(7) + v + B
    print('M =', M)
    print('scalar part:', M.scalar())
    print('grade-1 part:', M.get_grade(1))
    print('grade-2 part:', M.get_grade(2))

    # -- Symbolic example --
    print('\n--- Symbolic coordinates ---')
    theta = symbols('theta', real=True)
    u = cos(theta) * e1 + sin(theta) * e2
    print('u(theta) =', u)
    coefs = u.blade_coefs([e1, e2, e3])
    print('coefs:', coefs)

    # Evaluate at a specific angle
    vals = [c.subs(theta, pi / 4) for c in coefs]
    print('at theta=pi/4:', vals)
    print('numeric:', [float(v) for v in vals])


if __name__ == '__main__':
    main()
