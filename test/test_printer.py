import pytest
from sympy import Symbol

from galgebra.printer import GaPrinter, GaLatexPrinter
from galgebra.ga import Ga


def test_latex_flg():
    g3d, e_1, e_2, e_3 = Ga.build('e*1|2|3')
    t = Symbol('theta')

    mv = t + 0*e_1

    # it shouldn't matter whether the latex printer is enabled, if we call
    # the non-latex one we should get a non latex result.
    assert GaPrinter().doprint(mv) == 'theta'
    GaLatexPrinter.redirect()
    try:
        assert GaPrinter().doprint(mv) == 'theta'
    finally:
        GaLatexPrinter.restore()
