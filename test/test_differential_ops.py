from sympy import symbols, S
import pytest

from galgebra.ga import Ga
from galgebra.mv import Dop
from galgebra.dop import Sdop, Pdop


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
        assert (1234 * ex + ga.grad) * v == 1234 * ex * v + ga.grad * v != 0
        # check subtraction distributes
        assert (laplacian - ga.grad) * v == laplacian * v - ga.grad * v != 0
        assert (laplacian - 1234567) * v == laplacian * v - 1234567 * v != 0
        assert (1234 * ex - ga.grad) * v == 1234 * ex * v - ga.grad * v != 0
        # check unary subtraction distributes
        assert (-ga.grad) * v == -(ga.grad * v) != 0
        # check division is associative
        assert v * (ga.rgrad / 2) == (v * ga.rgrad) / 2 != 0
        # check multiplication is associative
        assert (ex * ga.grad) * v == ex * (ga.grad * v) != 0
        assert (20 * ga.grad) * v == 20 * (ga.grad * v) != 0
        assert v * (ga.rgrad * ex) == (v * ga.rgrad) * ex != 0
        assert v * (ga.rgrad * 20) == (v * ga.rgrad) * 20 != 0
        assert (laplacian * ga.grad) * v == laplacian * (ga.grad * v) != 0
        # check wedge is associative
        assert (ex ^ ga.grad) ^ v == ex ^ (ga.grad ^ v) != 0
        assert (20 ^ ga.grad) ^ v == 20 ^ (ga.grad ^ v) != 0
        assert v ^ (ga.rgrad ^ ex) == (v ^ ga.rgrad) ^ ex != 0
        assert v ^ (ga.rgrad ^ 20) == (v ^ ga.rgrad) ^ 20 != 0

    def test_empty_dop(self):
        """ Test that dop with zero terms is equivalent to multiplying by zero """
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        v = ga.mv('v', 'vector', f=True)

        make_zero = ga.dop([])
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

        # test comparison
        assert ga.grad == ga.grad
        assert not (ga.grad != ga.grad)
        assert ga.grad != laplacian
        assert not (ga.grad == laplacian)
        assert ga.grad != object()
        assert not (ga.grad == object())

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


class TestSdop(object):

    def test_deprecation(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        with pytest.warns(DeprecationWarning):
            ga.sPds

    def test_shorthand(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        assert Sdop(x) == Sdop([(S(1), Pdop({x: 1}))])

    def test_empty_sdop(self):
        """ Test that sdop with zero terms is equivalent to multiplying by zero """
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        v = ga.mv('v', 'vector', f=True)

        make_zero = Sdop([])
        assert make_zero * v == 0
        assert make_zero * make_zero * v == 0
        assert (make_zero + make_zero) * v == 0
        assert (-make_zero) * v == 0

    def test_associativity_and_distributivity(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        v = ga.mv('v', 'vector', f=True)
        laplacian = Sdop((ga.grad * ga.grad).terms)

        # check addition distributes
        assert (laplacian + 20) * v == laplacian * v + 20 * v != 0
        assert (20 + laplacian) * v == laplacian * v + 20 * v != 0
        # check subtraction distributes
        assert (laplacian - 20) * v == laplacian * v - 20 * v != 0
        assert (20 - laplacian) * v == 20 * v - laplacian * v != 0
        # check unary subtraction distributes
        assert (-laplacian) * v == -(laplacian * v) != 0
        # check multiplication is associative
        assert (20 * laplacian) * v == 20 * (laplacian * v) != 0
        assert (laplacian * laplacian) * v == laplacian * (laplacian * v) != 0


    def test_misc(self):
        """ Other miscellaneous tests """
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        laplacian = ga.sdop((ga.grad * ga.grad).terms)
        lap2 = laplacian*laplacian

        # test comparison
        assert lap2 == lap2
        assert not (lap2 != lap2)
        assert lap2 != laplacian
        assert not (lap2 == laplacian)
        assert lap2 != object()
        assert not (lap2 == object())


class TestPdop(object):

    def test_deprecation(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        # passing `None` is a deprecate way to spell `{}`
        with pytest.warns(DeprecationWarning):
            p = Pdop(None)
        assert p == Pdop({})

        with pytest.warns(DeprecationWarning):
            ga.Pdop_identity
        with pytest.warns(DeprecationWarning):
            ga.Pdiffs

    def test_misc(self):
        """ Other miscellaneous tests """
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        pxa = ga.pdop({x: 1})
        pxb = ga.pdop({x: 1})
        p1 = ga.pdop({})

        # test comparison
        assert pxa == pxb
        assert not (pxa != pxb)
        assert p1 != pxa
        assert not (p1 == pxa)
        assert pxa != object()
        assert not (pxa == object())
