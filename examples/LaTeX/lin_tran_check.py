from __future__ import print_function
from sympy import symbols, sin, cos, simplify
from galgebra.ga import Ga
from galgebra.printer import Format, xpdf, Eprint, Print_Function, Get_Program, latex
from galgebra.lt import Symbolic_Matrix


def main():
    # Print_Function()

    (x, y, z) = xyz = symbols('x,y,z',real=True)
    (o3d, ex, ey, ez) = Ga.build('e_x e_y e_z', g=[1, 1, 1], coords=xyz)
    grad = o3d.grad

    (u, v) = uv = symbols('u,v',real=True)
    (g2d, eu, ev) = Ga.build('e_u e_v', coords=uv)
    grad_uv = g2d.grad

    v_xyz = o3d.mv('v','vector')
    A_xyz = o3d.mv('A','vector',f=True)
    A_uv = g2d.mv('A','vector',f=True)

    print('#3d orthogonal ($A$ is vector function)')
    print('A =', A_xyz)
    print('%A^{2} =', A_xyz * A_xyz)
    print('grad|A =', grad | A_xyz)
    print('grad*A =', grad * A_xyz)

    print('v|(grad*A) =',v_xyz|(grad*A_xyz))

    print('#2d general ($A$ is vector function)')
    print('A =', A_uv)
    print('%A^{2} =', A_uv * A_uv)
    print('grad|A =', grad_uv | A_uv)
    print('grad*A =', grad_uv * A_uv)

    A = o3d.lt('A')

    print('#3d orthogonal ($A,\\;B$ are linear transformations)')
    print('A =', A)
    print(r'\f{mat}{A} =', A.matrix())
    print('\\f{\\det}{A} =', A.det())
    print('\\overline{A} =', A.adj())
    print('\\f{\\Tr}{A} =', A.tr())
    print('\\f{A}{e_x^e_y} =', A(ex^ey))
    print('\\f{A}{e_x}^\\f{A}{e_y} =', A(ex)^A(ey))

    B = o3d.lt('B')

    print('g =', o3d.g)
    print('%g^{-1} =', latex(o3d.g_inv))


    print('A + B =', A + B)
    print('AB =', A * B)
    print('A - B =', A - B)

    print('General Symmetric Linear Transformation')
    Asym = o3d.lt('A',mode='s')
    print('A =', Asym)
    print('General Antisymmetric Linear Transformation')
    Aasym = o3d.lt('A',mode='a')
    print('A =', Aasym)

    print('#2d general ($A,\\;B$ are linear transformations)')

    A2d = g2d.lt('A')

    print('g =', g2d.g)
    print('%g^{-1} =', latex(g2d.g_inv))
    print('%gg^{-1} =', latex(simplify(g2d.g * g2d.g_inv)))

    print('A =', A2d)
    print(r'\f{mat}{A} =', A2d.matrix())
    print('\\f{\\det}{A} =', A2d.det())
    A2d_adj = A2d.adj()
    print('\\overline{A} =', A2d_adj)
    print('\\f{mat}{\\overline{A}} =', latex(simplify(A2d_adj.matrix())))
    print('\\f{\\Tr}{A} =', A2d.tr())
    print('\\f{A}{e_u^e_v} =', A2d(eu^ev))
    print('\\f{A}{e_u}^\\f{A}{e_v} =', A2d(eu)^A2d(ev))

    B2d = g2d.lt('B')



    print('B =', B2d)
    print('A + B =', A2d + B2d)
    print('A - B =', A2d - B2d)

    # TODO: add this back when we drop Sympy 1.3. The 64kB of output is far too
    # printer-dependent
    if False:
        print('AB =', A2d * B2d)

    a = g2d.mv('a','vector')
    b = g2d.mv('b','vector')

    print(r'a|\f{\overline{A}}{b}-b|\f{\underline{A}}{a} =',((a|A2d.adj()(b))-(b|A2d(a))).simplify())

    m4d = Ga('e_t e_x e_y e_z', g=[1, -1, -1, -1],coords=symbols('t,x,y,z',real=True))

    T = m4d.lt('T')

    print('g =', m4d.g)

    print(r'\underline{T} =',T)
    print(r'\overline{T} =',T.adj())

    print(r'\f{\det}{\underline{T}} =',T.det())
    print(r'\f{\mbox{tr}}{\underline{T}} =',T.tr())

    a = m4d.mv('a','vector')
    b = m4d.mv('b','vector')

    print(r'a|\f{\overline{T}}{b}-b|\f{\underline{T}}{a} =',((a|T.adj()(b))-(b|T(a))).simplify())

    coords = (r, th, phi) = symbols('r,theta,phi', real=True)

    (sp3d, er, eth, ephi) = Ga.build('e_r e_th e_ph', g=[1, r**2, r**2*sin(th)**2], coords=coords)
    grad = sp3d.grad

    sm_coords = (u, v) = symbols('u,v', real=True)

    smap = [1, u, v]  # Coordinate map for sphere of r = 1

    sph2d = sp3d.sm(smap,sm_coords,norm=True)
    (eu, ev) = sph2d.mv()
    grad_uv = sph2d.grad

    F = sph2d.mv('F','vector',f=True)
    f = sph2d.mv('f','scalar',f=True)

    print('f =',f)
    print('grad*f =',grad_uv * f)

    print('F =',F)
    print('grad*F =',grad_uv * F)

    tp = (th,phi) = symbols('theta,phi',real=True)

    smap = [sin(th)*cos(phi),sin(th)*sin(phi),cos(th)]

    sph2dr = o3d.sm(smap,tp,norm=True)
    (eth, ephi) = sph2dr.mv()
    grad_tp = sph2dr.grad

    F = sph2dr.mv('F','vector',f=True)
    f = sph2dr.mv('f','scalar',f=True)

    print('f =',f)
    print('grad*f =',grad_tp * f)

    print('F =',F)
    print('grad*F =',grad_tp * F)

    return

def dummy():
    return

if __name__ == "__main__":
    #Eprint()
    Format()
    # Get_Program()
    main()
    # xpdf()
    xpdf(pdfprog=None)
