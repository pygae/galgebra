import unittest

from sympy import symbols
from galgebra.ga import Ga

class TestLt(unittest.TestCase):

    # reproduce gh-105
    def test_lt_matrix(self):
        base = Ga('a b', g=[1, 1], coords=symbols('x, y', real=True))
        a, b = base.mv()
        A = base.lt([a+b, 2*a-b])
        assert str(A) == 'Lt(a) = a + b\nLt(b) = 2*a - b'
        assert str(A.matrix()) == 'Matrix([[1, 2], [1, -1]])'
