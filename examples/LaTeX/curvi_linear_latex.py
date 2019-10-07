from __future__ import print_function
import sys

from sympy import symbols,sin,cos,sinh,cosh
from galgebra.ga import Ga
from galgebra.printer import Format, xpdf, Get_Program, Print_Function, Eprint

def derivatives_in_spherical_coordinates():
    #Print_Function()
    coords = (r,th,phi) = symbols('r theta phi', real=True)
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


def dummy():
    return


def main():
    #Get_Program()
    #Eprint()
    Format()
    derivatives_in_spherical_coordinates()
    derivatives_in_paraboloidal_coordinates()
    # FIXME This takes ~600 seconds
    # derivatives_in_elliptic_cylindrical_coordinates()
    derivatives_in_prolate_spheroidal_coordinates()
    #derivatives_in_oblate_spheroidal_coordinates()
    #derivatives_in_bipolar_coordinates()
    #derivatives_in_toroidal_coordinates()

    # xpdf()
    xpdf(pdfprog=None)
    return

if __name__ == "__main__":
    main()
