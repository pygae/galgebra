from __future__ import annotations

import os
import unittest
import pytest
import subprocess
from pathlib import Path

# Make SymPy available to this program:
import sympy 
from sympy import symbols

# Make GAlgebra available to this program:
from galgebra.ga import Ga  
# from galgebra.mv import *
from galgebra.printer import Fmt, GaPrinter, Format
# Fmt:       sets the way that a multivector's basis expansion is output.
# GaPrinter: makes GA output a little more readable.
# Format:    turns on latex printer.
from galgebra.gprinter import gFormat, gprint, gxpdf

ROOT = Path(__file__).parent.parent
DIR_TEST = ROOT / 'test'
DIR_GENERATED = DIR_TEST / 'generated'
DIR_FIXTURES = DIR_TEST / 'fixtures'
DIR_DIFF = DIR_TEST / 'diff'


class TestGprinter(unittest.TestCase):
    # @pytest.mark.skipif("TEST_GXPDF" not in os.environ, reason="Only run if TEST_GXPDF is set")
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

        # mkdir -p
        os.makedirs(DIR_GENERATED, exist_ok=True)
        os.makedirs(DIR_DIFF, exist_ok=True)
        os.chdir(DIR_GENERATED)

        gxpdf('test_gprinter.tex', pdfprog='tectonic', rm=False, null=False, evince=False)

        os.chdir(ROOT)
        retcode = subprocess.call([
            'diff-pdf',
            '--output-diff=%s' % (DIR_DIFF / 'test_gprinter-diff.pdf'),
            DIR_FIXTURES / 'test_gprinter.pdf',
            DIR_GENERATED / 'test_gprinter.pdf'
        ])
        assert retcode == 0

        # os.system('pdf2svg test_gprinter.pdf test_gprinter.svg')
