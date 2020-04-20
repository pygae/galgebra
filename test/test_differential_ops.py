import operator

from sympy import symbols, S
import pytest

from galgebra.ga import Ga
from galgebra.mv import Dop, Mv
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

    def test_constructor_errors(self):
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1])

        # list lengths must match
        with pytest.raises(ValueError, match='same length'):
            Dop([ex], [], ga=ga)

        # the two conventions can't be mixed
        mixed_args = [
            (ex, Pdop({})),
            (Sdop([]), ex),
        ]
        with pytest.raises(TypeError, match='pairs'):
            Dop(mixed_args, ga=ga)

        # ga must be non-none
        with pytest.raises(ValueError, match='must not be None'):
            Dop([], ga=None)

        # too few arguments
        with pytest.raises(TypeError, match='0 were given'):
            Dop(ga=ga)
        # too many arguments
        with pytest.raises(TypeError, match='3 were given'):
            Dop(1, 2, 3, ga=ga)


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

    def chain_with_pdop(self):
        x, y = symbols('x y', real=True)
        s = Sdop([(x, Pdop(x)), (y, Pdop(y))])
        # right-multiplication by Pdop chains only the pdiffs
        sp = s * Pdop(x)
        assert sp == Sdop([(x, Pdop(x) * Pdop(x)), (y, Pdop(y) * Pdop(x))])
        # left-multiplcation by Pdop invokes the product rule
        ps = Pdop(x) * s
        assert ps == Sdop([(x, Pdop(x) * Pdop(x)), (1, Pdop(x)), (y, Pdop(y) * Pdop(x))])

        # implicit multiplication
        assert ps == Pdop(x)(s)
        assert sp == s(Pdop(x))

        # no-op pdop
        assert s == Pdop({})(s)
        assert s == s(Pdop({}))

    def chain_with_mv(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)
        s = Sdop([(x, Pdop(x)), (y, Pdop(y))])

        assert type(ex * s) is Sdop
        assert type(s * ex) is Mv

        # type should be preserved even when the result is 0
        assert type(ex * Sdop([])) is Sdop
        assert type(Sdop([]) * ex) is Mv

        # As discussed with brombo, these operations are not well defined - if
        # you need them, you should be using `Dop` not `Sdop`.
        for op in [operator.xor, operator.or_, operator.lt, operator.gt]:
            with pytest.raises(TypeError):
                op(ex, s)
            with pytest.raises(TypeError):
                op(s, ex)

    def test_constructor_errors(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        # list lengths must match
        with pytest.raises(ValueError, match='same length'):
            Sdop([ex], [])

        # not a symbol or list
        with pytest.raises(TypeError, match='symbol or sequence is required'):
            Sdop(1)

        # not a pair of lists
        with pytest.raises(TypeError, match='must be lists'):
            Sdop([], 1)

        # too few arguments
        with pytest.raises(TypeError, match='0 were given'):
            Sdop()
        # too many arguments
        with pytest.raises(TypeError, match='3 were given'):
            Sdop(1, 2, 3)


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

    def test_multiply(self):
        coords = x, y, z = symbols('x y z', real=True)
        ga, ex, ey, ez = Ga.build('e*x|y|z', g=[1, 1, 1], coords=coords)

        p = Pdop(x)
        assert x * p == Sdop([(x, p)])
        assert ex * p == Sdop([(ex, p)])

        assert p * x == p(x) == S(1)
        assert p * ex == p(ex) == S(0)
        assert type(p(ex)) is Mv

        # These are not defined for consistency with Sdop
        for op in [operator.xor, operator.or_, operator.lt, operator.gt]:
            with pytest.raises(TypeError):
                op(ex, p)
            with pytest.raises(TypeError):
                op(p, ex)

    def test_constructor_errors(self):
        # not a symbol or dict
        with pytest.raises(TypeError, match='dictionary or symbol is required'):
            Pdop(1)
