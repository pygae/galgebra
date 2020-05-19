import contextlib
import io
import textwrap
import sys

import pytest
from sympy import Symbol

from galgebra.printer import GaPrinter, GaLatexPrinter, oprint
from galgebra.ga import Ga


has_ordered_dictionaries = sys.version_info >= (3, 6)


def test_latex_flg_GaPrinter():
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


def test_latex_flg_Symbol_sortkey():
    """
    Symbol.sort_key should not be affected by whether latex printing is enabled.

    `sort_key` affects the order of expressions returned by `sympy.simplify`.
    """
    t = Symbol('theta')

    t_sort = t.sort_key()

    # sort key is cached, and we want to compare the behavior between the cache
    # being populated from within a latex context to it being populated from
    # outside it
    t.sort_key.cache_clear()

    GaLatexPrinter.redirect()
    try:
        t_latex_sort = t.sort_key()
    finally:
        GaLatexPrinter.restore()

    assert t_sort == t_latex_sort


def test_oprint():
    s = io.StringIO()
    with contextlib.redirect_stdout(s):
        oprint(
            'int', 1,
            'dictionary', dict(a=1, b=2),
            'set', {1},
            'tuple', (1, 2),
            'list', [1, 2, 3],
            'str', 'a quote: "',
        )

    if has_ordered_dictionaries:
        assert s.getvalue() == textwrap.dedent("""\
            int        = 1
            dictionary = {a:1,b:2}
            set        = {1}
            tuple      = (1,2)
            list       = [1,2,3]
            str        = a quote: "
            """)


def test_oprint_dict_mode():
    s = io.StringIO()
    with contextlib.redirect_stdout(s):
        oprint(
            'int', 1,
            'dictionary', dict(a=1, b=2),
            'set', {1},
            'tuple', (1, 2),
            'list', [1, 2, 3],
            'str', 'a quote: "',
            dict_mode=True
        )

    if has_ordered_dictionaries:
        assert s.getvalue() == textwrap.dedent("""\
            int   = 1
            dictionary:
            a -> 1
            b -> 2
            set   = {1}
            tuple = (1,2)
            list  = [1,2,3]
            str   = a quote: "
            """)
