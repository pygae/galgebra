import unittest

import pytest
from sympy import symbols, S, Matrix, Symbol

from galgebra.ga import Ga
from galgebra.lt import Mlt


class TestLt(unittest.TestCase):

    # reproduce gh-540: callable Lt must accept zero (e.g. projection maps)
    def test_lt_callable_zero(self):
        ga, e1, e2, e3 = Ga.build('e*1|2|3', g=[1, 1, 1])
        # projection onto e1: maps e2 and e3 to zero
        L = ga.lt(lambda x: (x | e1) * e1)
        assert L(e1) == e1
        assert L(e2).is_zero()
        assert L(e3).is_zero()
        # zero map: every basis vector maps to zero
        L_zero = ga.lt(lambda x: x - x)
        for basis_v in [e1, e2, e3]:
            assert L_zero(basis_v).is_zero()

    # reproduce gh-105
    def test_lt_matrix(self):
        base = Ga('a b', g=[1, 1], coords=symbols('x, y', real=True))
        a, b = base.mv()
        A = base.lt([a+b, 2*a-b])
        assert str(A) == 'Lt(a) = a + b\nLt(b) = 2*a - b'
        assert str(A.matrix()) == 'Matrix([[1, 2], [1, -1]])'

    # reproduce gh-461: lt.matrix() on non-Euclidean metrics
    def test_lt_matrix_oblique(self):
        # oblique metric g=[[1,1],[1,2]]: matrix() must not include metric factors
        coords = symbols('x y', real=True)
        ga, e1, e2 = Ga.build('e*1|2', g=[[1, 1], [1, 2]], coords=coords)
        L = ga.lt([e1 + e2, 2*e1 - e2])
        assert L.matrix() == Matrix([[1, 2], [1, -1]])

        # Minkowski metric g=diag(1,-1): same check
        ga2, f0, f1 = Ga.build('e*0|1', g=[1, -1], coords=symbols('t x', real=True))
        L2 = ga2.lt([f0 + f1, 2*f0 - f1])
        assert L2.matrix() == Matrix([[1, 2], [1, -1]])

    # reproduce gh-461: Ga.lt('f') on oblique metric must not mix in metric tensor
    def test_lt_generic_oblique(self):
        coords = symbols('x y', real=True)
        ga, e1, e2 = Ga.build('e*1|2', g=[[1, 1], [1, 2]], coords=coords)
        F = ga.lt('f')
        c1 = F(e1).get_coefs(1)
        c2 = F(e2).get_coefs(1)
        # each coefficient must be a plain symbol, not a metric-weighted expression
        assert all(isinstance(s, Symbol) for s in c1 + c2)
        # and the matrix must match those coefficients exactly
        assert F.matrix() == Matrix([[c1[0], c2[0]], [c1[1], c2[1]]])

    def test_lt_function(self):
        """ Test construction from a function """
        base = Ga('a b', g=[1, 1], coords=symbols('x, y', real=True))
        a, b = base.mv()

        def not_linear(x):
            return x * x
        with pytest.raises(ValueError, match='linear'):
            base.lt(not_linear)

        def not_vector(x):
            return x + S.One
        with pytest.raises(ValueError, match='vector'):
            base.lt(not_vector)

        def ok(x):
            return (x | b) * a + 2*x
        f = base.lt(ok)
        x = base.mv('x', 'vector')
        y = base.mv('y', 'vector')
        assert f(x) == ok(x)
        assert f(x^y) == ok(x)^ok(y)
        assert f(1 + 2*(x^y)) == 1 + 2*(ok(x)^ok(y))

    def test_deprecations(self):
        base = Ga('a b', g=[1, 1], coords=symbols('x, y', real=True))
        l = base.lt([[1, 2], [3, 4]])
        with pytest.warns(DeprecationWarning):
            assert l.X == l.Ga.coord_vec
        with pytest.warns(DeprecationWarning):
            assert l.coords == l.Ga.coords

        l = base.lt('L', mode='a')
        with pytest.warns(DeprecationWarning):
            assert l.mode == 'a'
        with pytest.warns(DeprecationWarning):
            assert not l.fct_flg

        l = base.lt('L', mode='s', f=True)
        with pytest.warns(DeprecationWarning):
            assert l.mode == 's'
        with pytest.warns(DeprecationWarning):
            assert l.fct_flg



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

    def test_str_arithmetic(self):
        """Mlt arithmetic on string-constructed tensors routes through the
        component-expression constructor and must not raise NotImplementedError."""
        coords = symbols('x y', real=True)
        g, e1, e2 = Ga.build('e*1|2', coords=coords, g=[1, 1])

        a1 = g.mv('a1', 'vector')
        a2 = g.mv('a2', 'vector')

        S = Mlt('S', g, nargs=2)
        T = Mlt('T', g, nargs=2)

        assert (S + T)(a1, a2) == S(a1, a2) + T(a1, a2)
        assert (S - T)(a1, a2) == S(a1, a2) - T(a1, a2)

    def test_from_component_expression(self):
        """Mlt constructed from a pre-built component expression (issue #578)."""
        coords = symbols('x y', real=True)
        g, e1, e2 = Ga.build('e*1|2', coords=coords, g=[1, 1])

        a1 = g.mv('a1', 'vector')
        a2 = g.mv('a2', 'vector')

        # Build the same rank-2 tensor as test_from_str but pass the
        # fvalue expression directly instead of a string name.
        T_str = Mlt('T', g, nargs=2)
        fvalue = T_str.fvalue  # sympy expression in slot variables

        # Reconstruct using the component-expression path
        T_expr = Mlt(fvalue, g, nargs=2)
        assert T_expr(a1, a2) == T_str(a1, a2)

        # nargs is required for component expressions
        with pytest.raises(TypeError):
            Mlt(fvalue, g)

    def test_deprecations(self):
        g = Ga('e*a|b', g=[1, 1])
        with pytest.warns(DeprecationWarning):
            assert Mlt.extact_basis_indexes(g) == ['a', 'b']
