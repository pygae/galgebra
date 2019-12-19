from sympy import symbols
import pytest

from galgebra.ga import Ga
from galgebra.mv import Dop


class TestDop(object):

    def test_associativity_and_distributivity(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        v = ga.mv('v', 'vector', f=True)
        laplacian = ga.grad * ga.grad
        rlaplacian = ga.rgrad * ga.rgrad

        # check addition distributes
        assert (laplacian + ga.grad) * v == laplacian * v + ga.grad * v != 0
        assert (laplacian + 1234567) * v == laplacian * v + 1234567 * v != 0
        # check subtraction distributes
        assert (laplacian - ga.grad) * v == laplacian * v - ga.grad * v != 0
        assert (laplacian - 1234567) * v == laplacian * v - 1234567 * v != 0
        # check unary subtraction distributes
        assert (-ga.grad) * v == -(ga.grad * v) != 0
        # check division is associative
        assert v * (ga.rgrad / 2) == (v * ga.rgrad) / 2 != 0
        # check multiplication is associative
        assert (20 * ga.grad) * v == 20 * (ga.grad * v) != 0
        assert v * (ga.rgrad * 20) == (v * ga.rgrad) * 20 != 0
        assert (ex ^ ga.grad) ^ v == ex ^ (ga.grad ^ v) != 0
        assert (20 ^ ga.grad) ^ v == 20 ^ (ga.grad ^ v) != 0
        assert v ^ (ga.rgrad ^ ex) == (v ^ ga.rgrad) ^ ex != 0
        assert v ^ (ga.rgrad ^ 20) == (v ^ ga.rgrad) ^ 20 != 0

    def test_empty_dop(self):
        """ Test that dop with zero terms is equivalent to multiplying by zero """
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        v = ga.mv('v', 'vector', f=True)

        make_zero = Dop([], ga=ga)
        assert make_zero * v == 0
        assert make_zero * make_zero * v == 0
        assert (make_zero + make_zero) * v == 0
        assert (-make_zero) * v == 0

    def test_misc(self):
        """ Other miscellaneous tests """
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        v = ga.mv('v', 'vector', f=True)
        laplacian = ga.grad * ga.grad
        rlaplacian = ga.rgrad * ga.rgrad

        # laplacian is a scalar operator, so applying it from either side
        # is the same
        assert laplacian * v == v * rlaplacian

        assert laplacian.is_scalar()
        assert not ga.grad.is_scalar()

        assert ga.grad == ga.grad
        assert ga.grad != laplacian

        # inconsistent cmpflg, not clear which side the operator goes on
        with pytest.raises(ValueError):
            ga.grad + ga.rgrad
        with pytest.raises(ValueError):
            ga.grad * ga.rgrad

    def test_mixed_algebras(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga1, ex1, ey1, ez1 = Ga.build('e1*x|y|z', g=[1, 1, 1], coords=coords)
        ga2, ex2, ey2, ez2 = Ga.build('e2*x|y|z', g=[1, 1, 1], coords=coords)
        assert ga1 != ga2

        v1 = ga1.mv('v', 'vector', f=True)
        v2 = ga2.mv('v', 'vector', f=True)

        with pytest.raises(ValueError):
            ga1.grad * v2
        with pytest.raises(ValueError):
            v1 * ga2.rgrad
        with pytest.raises(ValueError):
            ga1.grad * ga2.grad

    def test_components(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        components = ga.grad.components()
        assert components == (
            ex * (ex | ga.grad),
            ey * (ey | ga.grad),
            ez * (ez | ga.grad),
        )
