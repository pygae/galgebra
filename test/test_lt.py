import unittest

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

        # calling the Mlt is like calling the function
        assert T(a1, a2) == TA(a1, a2)
