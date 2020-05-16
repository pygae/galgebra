import unittest

import pytest
from sympy import symbols

from galgebra.ga import Ga
from galgebra.lt import Mlt


class TestLt(unittest.TestCase):

    # reproduce gh-105
    def test_lt_matrix(self):
        base = Ga('a b', g=[1, 1], coords=symbols('x, y', real=True))
        a, b = base.mv()
        A = base.lt([a+b, 2*a-b])
        assert str(A) == 'Lt(a) = a + b\nLt(b) = 2*a - b'
        assert str(A.matrix()) == 'Matrix([[1, 2], [1, -1]])'


class TestMlt(unittest.TestCase):

    def test_basic(self):
        # from TensorDef.py
        coords = symbols('t x y z', real=True)
        st4d, g0, g1, g2, g3 = Ga.build('gamma*t|x|y|z', g=[1, -1, -1, -1],
                                        coords=coords)

        A = st4d.mv('T', 'bivector')

        def TA(a1, a2):
            return A | (a1 ^ a2)

        T = Mlt(TA, st4d)

        # tests begin

        a1 = st4d.mv('a1', 'vector')
        a2 = st4d.mv('a2', 'vector')
        a3 = st4d.mv('a3', 'vector')
        a4 = st4d.mv('a4', 'vector')

        # calling the Mlt is like calling the function
        assert T(a1, a2) == TA(a1, a2)

        # for addition, argument slots are reused
        assert (T + T)(a1, a2) == T(a1, a2) + T(a1, a2)
        assert (T - T)(a1, a2) == T(a1, a2) - T(a1, a2)

        # for multiplication, argument slots are chained
        assert (T * T)(a1, a2, a3, a4) == TA(a1, a2) * T(a3, a4)
        assert (T ^ T)(a1, a2, a3, a4) == TA(a1, a2) ^ T(a3, a4)
        assert (T | T)(a1, a2, a3, a4) == TA(a1, a2) | T(a3, a4)

        # Test linearity properties. Note that this behavior is implied by our
        # test that T and TA are equivalent above, but it does exercise
        # `Mlt.__call__` with compound expressions as arguments.
        alpha = st4d.mv('alpha', 'scalar')
        b = st4d.mv('b', 'vector')

        assert T(alpha * a1, a2) == alpha * T(a1, a2)
        assert T(a1, alpha * a2) == alpha * T(a1, a2)
        assert T(a1 + b, a2) == T(a1, a2) + T(b, a2)
        assert T(a1, a2 + b) == T(a1, a2) + T(a1, b)

    def test_from_str(self):
        coords = symbols('x y', real=True)
        g, e1, e2 = Ga.build('e*1|2', coords=coords, g=[1, 1])

        a1 = g.mv('a1', 'vector')
        a2 = g.mv('a2', 'vector')
        a1x, a1y = a1.get_coefs(1)
        a2x, a2y = a2.get_coefs(1)

        # one-d
        T = Mlt('T', g, nargs=1)
        v = T(a1)

        # Two new symbols created
        Tx, Ty = sorted(v.free_symbols - {a1x, a1y}, key=lambda x: x.sort_key())
        assert v == (
            Tx * a1x +
            Ty * a1y
        )

        # two-d
        T = Mlt('T', g, nargs=2)
        v = T(a1, a2)

        # four new symbols created
        Txx, Txy, Tyx, Tyy = sorted(v.free_symbols - {a1x, a1y, a2x, a2y}, key=lambda x: x.sort_key())
        assert v == (
            Txx * a1x * a2x +
            Txy * a1x * a2y +
            Tyx * a1y * a2x +
            Tyy * a1y * a2y
        )

    def test_deprecations(self):
        g = Ga('e*a|b', g=[1, 1])
        with pytest.warns(DeprecationWarning):
            assert Mlt.extact_basis_indexes(g) == ['a', 'b']
