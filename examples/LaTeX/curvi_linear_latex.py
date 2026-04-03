from __future__ import print_function
import sys
from contextlib import contextmanager

from sympy import symbols,sin,cos,sinh,cosh,trigsimp
from galgebra.ga import Ga
from galgebra.metric import Simp
from galgebra.printer import Format, xpdf, Print_Function, Eprint


@contextmanager
def _no_simp_build():
    """Use an identity simplifier during Ga.build, then restore the caller's Simp profile.

    ``Ga.build(norm=True)`` calls ``Simp.apply`` ~70 times during metric
    normalisation.  Since SymPy 1.13 (PR #26390) each call traverses every
    trig/hyperbolic node via ``.replace(Mul, ...)`` inside ``TR3``/``futrig``,
    which is always a no-op for symbolic trig arguments but costs ~25 s per call.

    By substituting an identity simplifier for the build phase we skip all ~70
    traversals.  The callers' ``Simp.modes`` (set by ``main()`` to
    ``trigsimp(method='old')``) remain active for the output expressions, so
    ``_sympystr`` still simplifies each printed result once when formatting.
    """
    orig = Simp.modes[:]
    Simp.profile([lambda e: e])  # identity — no simplification during build
    try:
        yield
    finally:
        Simp.profile(orig)


def derivatives_in_spherical_coordinates():
    #Print_Function()
    coords = (r,th,phi) = symbols('r theta phi', real=True)
    with _no_simp_build():
        (sp3d,er,eth,ephi) = Ga.build('e_r e_theta e_phi',g=[1,r**2,r**2*sin(th)**2],coords=coords)
    grad = sp3d.grad

    f = sp3d.mv('f','scalar',f=True)
    A = sp3d.mv('A','vector',f=True)
    B = sp3d.mv('B','bivector',f=True)

    print('#Derivatives in Spherical Coordinates')

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    print('grad\\times A = -I*(grad^A) =',-sp3d.i*(grad^A))
    print('%\\nabla^{2}f =',grad|(grad*f))
    print('grad^B =',grad^B)

    """
    print '( \\nabla\\W\\nabla )\\bm{e}_{r} =',((grad^grad)*er).trigsimp()
    print '( \\nabla\\W\\nabla )\\bm{e}_{\\theta} =',((grad^grad)*eth).trigsimp()
    print '( \\nabla\\W\\nabla )\\bm{e}_{\\phi} =',((grad^grad)*ephi).trigsimp()
    """

    return


def derivatives_in_paraboloidal_coordinates():
    #Print_Function()
    coords = (u,v,phi) = symbols('u v phi', real=True)
    with _no_simp_build():
        (par3d,er,eth,ephi) = Ga.build('e_u e_v e_phi',X=[u*v*cos(phi),u*v*sin(phi),(u**2-v**2)/2],coords=coords,norm=True)
    grad = par3d.grad

    f = par3d.mv('f','scalar',f=True)
    A = par3d.mv('A','vector',f=True)
    B = par3d.mv('B','bivector',f=True)

    print('#Derivatives in Paraboloidal Coordinates')

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    (-par3d.i*(grad^A)).Fmt(3,'grad\\times A = -I*(grad^A)')
    print('grad^B =',grad^B)

    return


def derivatives_in_elliptic_cylindrical_coordinates():
    #Print_Function()
    a = symbols('a', real=True)
    coords = (u,v,z) = symbols('u v z', real=True)
    with _no_simp_build():
        (elip3d,er,eth,ephi) = Ga.build('e_u e_v e_z',X=[a*cosh(u)*cos(v),a*sinh(u)*sin(v),z],coords=coords,norm=True)
    grad = elip3d.grad

    f = elip3d.mv('f','scalar',f=True)
    A = elip3d.mv('A','vector',f=True)
    B = elip3d.mv('B','bivector',f=True)

    print('#Derivatives in Elliptic Cylindrical Coordinates')

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    print('-I*(grad^A) =',-elip3d.i*(grad^A))
    print('grad^B =',grad^B)
    return


def derivatives_in_prolate_spheroidal_coordinates():
    #Print_Function()
    a = symbols('a', real=True)
    coords = (xi,eta,phi) = symbols('xi eta phi', real=True)
    with _no_simp_build():
        (ps3d,er,eth,ephi) = Ga.build('e_xi e_eta e_phi',X=[a*sinh(xi)*sin(eta)*cos(phi),a*sinh(xi)*sin(eta)*sin(phi),
                                                            a*cosh(xi)*cos(eta)],coords=coords,norm=True)
    grad = ps3d.grad

    f = ps3d.mv('f','scalar',f=True)
    A = ps3d.mv('A','vector',f=True)
    B = ps3d.mv('B','bivector',f=True)

    print('#Derivatives in Prolate Spheroidal Coordinates')

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    (-ps3d.i*(grad^A)).Fmt(3,'-I*(grad^A)')
    (grad^B).Fmt(3,'grad^B')
    return


def derivatives_in_oblate_spheroidal_coordinates():
    Print_Function()
    a = symbols('a', real=True)
    coords = (xi,eta,phi) = symbols('xi eta phi', real=True)
    with _no_simp_build():
        (os3d,er,eth,ephi) = Ga.build('e_xi e_eta e_phi',X=[a*cosh(xi)*cos(eta)*cos(phi),a*cosh(xi)*cos(eta)*sin(phi),
                                                            a*sinh(xi)*sin(eta)],coords=coords,norm=True)
    grad = os3d.grad

    f = os3d.mv('f','scalar',f=True)
    A = os3d.mv('A','vector',f=True)
    B = os3d.mv('B','bivector',f=True)

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    print('-I*(grad^A) =',-os3d.i*(grad^A))
    print('grad^B =',grad^B)
    return


def derivatives_in_bipolar_coordinates():
    Print_Function()
    a = symbols('a', real=True)
    coords = (u,v,z) = symbols('u v z', real=True)
    with _no_simp_build():
        (bp3d,eu,ev,ez) = Ga.build('e_u e_v e_z',X=[a*sinh(v)/(cosh(v)-cos(u)),a*sin(u)/(cosh(v)-cos(u)),z],coords=coords,norm=True)
    grad = bp3d.grad

    f = bp3d.mv('f','scalar',f=True)
    A = bp3d.mv('A','vector',f=True)
    B = bp3d.mv('B','bivector',f=True)

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    print('-I*(grad^A) =',-bp3d.i*(grad^A))
    print('grad^B =',grad^B)
    return


def derivatives_in_toroidal_coordinates():
    Print_Function()
    a = symbols('a', real=True)
    coords = (u,v,phi) = symbols('u v phi', real=True)
    with _no_simp_build():
        (t3d,eu,ev,ephi) = Ga.build('e_u e_v e_phi',X=[a*sinh(v)*cos(phi)/(cosh(v)-cos(u)),
                                                        a*sinh(v)*sin(phi)/(cosh(v)-cos(u)),
                                                        a*sin(u)/(cosh(v)-cos(u))],coords=coords,norm=True)
    grad = t3d.grad

    f = t3d.mv('f','scalar',f=True)
    A = t3d.mv('A','vector',f=True)
    B = t3d.mv('B','bivector',f=True)

    print('f =',f)
    print('A =',A)
    print('B =',B)

    print('grad*f =',grad*f)
    print('grad|A =',grad|A)
    print('-I*(grad^A) =',-t3d.i*(grad^A))
    print('grad^B =',grad^B)
    return


def main():
    #Eprint()
    Format()

    # SymPy >= 1.13 (PR #26390) added a slow O(N*M) traversal inside
    # sympy.simplify.fu that causes timeouts on curvilinear coordinate
    # expressions.  Use trigsimp(method='old') via Simp.profile to avoid
    # that code path entirely for this example.
    #
    # _no_simp_build() inside each function further avoids ~70 redundant
    # Simp.apply calls during Ga.build by substituting an identity simplifier
    # for the build phase only.  The trigsimp(method='old') profile below
    # then applies once per output expression via _sympystr at print time.
    orig_modes = Simp.modes[:]
    Simp.profile([lambda e: trigsimp(e, method='old')])
    try:
        derivatives_in_spherical_coordinates()
        derivatives_in_paraboloidal_coordinates()
        # FIXME This takes ~600 seconds
        # derivatives_in_elliptic_cylindrical_coordinates()
        derivatives_in_prolate_spheroidal_coordinates()
        #derivatives_in_oblate_spheroidal_coordinates()
        #derivatives_in_bipolar_coordinates()
        #derivatives_in_toroidal_coordinates()
    finally:
        Simp.profile(orig_modes)

    # xpdf()
    xpdf(pdfprog=None)
    return

if __name__ == "__main__":
    main()
