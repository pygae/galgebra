import unittest

import pytest
from sympy import Symbol, symbols
from galgebra.ga import Ga
from galgebra import utils

class TestLt(unittest.TestCase):

    # reproduce gh-105
    def test_lt_matrix(self):
        base = Ga('a b', g=[1,1], coords=symbols('x,y',real=True))
        a,b = base.mv()
        A = base.lt([a+b,a-b])
        assert str(A) == 'Lt(a) = a + b\nLt(b) = a - b'
        assert str(A.matrix()) == 'Matrix([[1, 1], [1, -1]])'