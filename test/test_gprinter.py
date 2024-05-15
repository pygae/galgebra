from __future__ import annotations

import os
import unittest
import pytest
import subprocess
from pathlib import Path

# Make SymPy available to this program:
import sympy 
from sympy import symbols, simplify, sin, cos

# Make GAlgebra available to this program:
from galgebra.ga import Ga  
# from galgebra.mv import *
from galgebra.printer import Fmt, GaPrinter, Format
# Fmt:       sets the way that a multivector's basis expansion is output.
# GaPrinter: makes GA output a little more readable.
# Format:    turns on latex gprinter.
from galgebra.gprinter import gFormat, gprint, gxpdf, gPrint_Function

ROOT = Path(__file__).parent.parent
DIR_TEST = ROOT / 'test'
DIR_GENERATED = DIR_TEST / 'generated'
DIR_FIXTURES = DIR_TEST / 'fixtures'
DIR_DIFF = DIR_TEST / 'diff'


def lin_tran():
    gFormat()
    gPrint_Function()
    (x, y, z) = xyz = symbols('x,y,z',real=True)
    (o3d, ex, ey, ez) = Ga.build('e_x e_y e_z', g=[1, 1, 1], coords=xyz)
    grad = o3d.grad
    (u, v) = uv = symbols('u,v',real=True)
    (g2d, eu, ev) = Ga.build('e_u e_v', coords=uv)
    grad_uv = g2d.grad
    v_xyz = o3d.mv('v','vector')
    A_xyz = o3d.mv('A','vector',f=True)
    A_uv = g2d.mv('A','vector',f=True)
    # gprint(r'\text{3d orthogonal ($A$ is vector function)}')
    gprint('A =', A_xyz)
    gprint('A^{2} =', A_xyz * A_xyz)
    gprint(r'\nabla \cdot A =', grad | A_xyz)
    gprint(r'\nabla A =', grad * A_xyz)
    gprint('v|(grad*A) =',v_xyz|(grad*A_xyz))
    gprint(r'\text{2d general ($A$ is vector function)}')
    gprint('A =', A_uv)
    gprint('A^{2} =', A_uv * A_uv)
    gprint(r'\nabla \cdot A =', grad_uv | A_uv)
    gprint(r'\nabla A =', grad_uv * A_uv)
    A = o3d.lt('A')
    gprint(r'\text{3d orthogonal ($A$, $B$ are linear transformations)}')
    gprint('A =', A)
    gprint(r'\f{mat}{A} =', A.matrix())
    gprint('\\f{\\det}{A} =', A.det())
    gprint('\\overline{A} =', A.adj())
    gprint('\\f{\\Tr}{A} =', A.tr())
    gprint('\\f{A}{e_x \\wedge e_y} =', A(ex^ey))
    gprint('\\f{A}{e_x} \\wedge \\f{A}{e_y} =', A(ex)^A(ey))
    B = o3d.lt('B')
    gprint('g =', o3d.g)
    gprint('g^{-1} =', o3d.g_inv)

    gprint('A + B =', A + B)
    gprint('AB =', A * B)
    gprint('A - B =', A - B)
    gprint('General Symmetric Linear Transformation')
    Asym = o3d.lt('A',mode='s')
    gprint('A =', Asym)
    gprint('General Antisymmetric Linear Transformation')
    Aasym = o3d.lt('A',mode='a')
    gprint('A =', Aasym)
    gprint(r'\text{2d general ($A,\\;B$ are linear transformations)}')
    A2d = g2d.lt('A')
    gprint('g =', g2d.g)
    gprint('g^{-1} =', g2d.g_inv)
    gprint('gg^{-1} =', simplify(g2d.g * g2d.g_inv))
    gprint('A =', A2d)
    gprint(r'\f{mat}{A} =', A2d.matrix())
    gprint('\\f{\\det}{A} =', A2d.det())
    A2d_adj = A2d.adj()
    gprint('\\overline{A} =', A2d_adj)
    gprint('\\f{mat}{\\overline{A}} =', simplify(A2d_adj.matrix()))
    gprint('\\f{\\Tr}{A} =', A2d.tr())
    gprint('\\f{A}{e_u \\wedge e_v} =', A2d(eu^ev))
    gprint('\\f{A}{e_u} \\wedge \\f{A}{e_v} =', A2d(eu)^A2d(ev))
    B2d = g2d.lt('B')

    gprint('B =', B2d)
    gprint('A + B =', A2d + B2d)
    gprint('A - B =', A2d - B2d)
    gprint('AB =', A2d * B2d)
    a = g2d.mv('a','vector')
    b = g2d.mv('b','vector')
    gprint(r'a|\f{\overline{A}}{b}-b|\f{\underline{A}}{a} =',((a|A2d.adj()(b))-(b|A2d(a))).simplify())
    m4d = Ga('e_t e_x e_y e_z', g=[1, -1, -1, -1],coords=symbols('t,x,y,z',real=True))
    T = m4d.lt('T')
    gprint('g =', m4d.g)
    gprint(r'\underline{T} =',T)
    gprint(r'\overline{T} =',T.adj())
    gprint(r'\f{\det}{\underline{T}} =',T.det())
    gprint(r'\f{\mbox{tr}}{\underline{T}} =',T.tr())
    a = m4d.mv('a','vector')
    b = m4d.mv('b','vector')
    gprint(r'a|\f{\overline{T}}{b}-b|\f{\underline{T}}{a} =',((a|T.adj()(b))-(b|T(a))).simplify())
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
    gprint('f =',f)
    gprint(r'\nabla f =',grad_uv * f)
    gprint('F =',F)
    gprint(r'\nabla F =',grad_uv * F)
    tp = (th,phi) = symbols('theta,phi',real=True)
    smap = [sin(th)*cos(phi),sin(th)*sin(phi),cos(th)]
    sph2dr = o3d.sm(smap,tp,norm=True)
    (eth, ephi) = sph2dr.mv()
    grad_tp = sph2dr.grad
    F = sph2dr.mv('F','vector',f=True)
    f = sph2dr.mv('f','scalar',f=True)
    gprint('f =',f)
    gprint(r'\nabla f =',grad_tp * f)
    gprint('F =',F)
    gprint(r'\nabla F =',grad_tp * F)
    return


class TestGprinter(unittest.TestCase):
    def setUp(self):
        os.makedirs(DIR_GENERATED, exist_ok=True)
        os.makedirs(DIR_DIFF, exist_ok=True)
        os.chdir(DIR_GENERATED)

    def gen_pdf(self, name, documentclass):
        gxpdf('%s.tex' % name, pdfprog='tectonic', rm=False, null=False, evince=False, documentclass=documentclass)

    def check_pdf(self, name, expected_retcode=0):
        retcode = subprocess.call([
            'diff-pdf',
            '--output-diff=%s' % (DIR_DIFF / ('%s-diff.pdf' % name)),
            DIR_FIXTURES / ('%s.pdf' % name),
            DIR_GENERATED / ('%s.pdf' % name)
        ])
        assert retcode == expected_retcode

    @pytest.mark.skipif("TEST_GXPDF" not in os.environ, reason="Only run if TEST_GXPDF is set")
    def test_gxpdf(self):
        gFormat()
        # Set up standard G^3 geometric algebra
        g3coords = (x, y, z) = symbols('x y z', real=True)  # Without real=True, symbols are complex
        g3 = Ga(r'\mathbf{e}', g=[1, 1, 1], coords=g3coords)
        (ex, ey, ez) = g3.mv()      # Program names of basis vectors.
        (exr, eyr, ezr) = g3.mvr()  # Program names of reciprocal basis vectors.

        B = g3.mv('B', 'bivector')
        Fmt(1)  # Set Fmt globally
        gprint(r'\mathbf{B} =', B)         # B will be bold.
        gprint(r'\mathbf{B} =', B.Fmt(3))  # Fmt(3) here only.
        gprint(r'\mathbf{B} =', B)         # Global Fmt remembered.

        gprint(r'\mathbf{B}^2 =', B*B)

        M = g3.mv('M', 'mv')
        gprint(r'\langle \mathbf{M} \rangle_2 =', M.grade(2)) 
        # grade(2) could be replaced by, e.g., odd(), or omitted altogether.

        gprint(r'\alpha_1\mathbf{X}/\gamma_r^3')

        grad = g3.grad

        gprint(r'{\nabla} = ', grad)

        self.gen_pdf('test_gprinter', documentclass='book')
        self.check_pdf('test_gprinter')

    # This test is a POC, and should skip for now
    @pytest.mark.skip
    def test_lt_pdf(self):
        lin_tran()
        self.gen_pdf('test_lt_pdf', documentclass='report')
        self.check_pdf('test_lt_pdf', expected_retcode=1)