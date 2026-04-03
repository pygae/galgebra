from packaging.version import Version

import pytest
import sympy
from sympy import symbols, S
from galgebra.ga import Ga
from galgebra.mv import proj, undual, g_invol, exp, norm, norm2, mag, mag2, ccon, rev, scalar, qform, sp, inv, shirokov_inverse, hitzer_inverse


class TestMv:

    def test_deprecations(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3')
        with pytest.warns(DeprecationWarning):
            ga.mv_x
        with pytest.warns(DeprecationWarning):
            ga.mv_I

    def test_is_base(self):
        """
        Various tests on several multivectors.
        """
        _g3d, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        assert (e_1).is_base()
        assert (e_2).is_base()
        assert (e_3).is_base()
        assert (e_1 ^ e_2).is_base()
        assert (e_2 ^ e_3).is_base()
        assert (e_1 ^ e_3).is_base()
        assert (e_1 ^ e_2 ^ e_3).is_base()

        assert not (2 * e_1).is_base()
        assert not (e_1 + e_2).is_base()
        assert not (e_3 * 4).is_base()
        assert not ((3 * e_1) ^ e_2).is_base()
        assert not (2 * (e_2 ^ e_3)).is_base()
        assert not (e_3 ^ e_1).is_base()
        assert not (e_2 ^ e_1 ^ e_3).is_base()

    def test_get_coefs(self):
        _g3d, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        assert (e_1 * 3 + e_3).get_coefs(1) == [3, 0, 1]

        # can always get coefficients of 0
        assert (0*e_1).get_coefs(0) == [0]
        assert (0*e_1).get_coefs(1) == [0, 0, 0]
        assert (0*e_1).get_coefs(2) == [0, 0, 0]
        assert (0*e_1).get_coefs(3) == [0]

        # grade is wrong
        with pytest.raises(ValueError):
            (e_1 ^ e_2).get_coefs(1)
        with pytest.raises(ValueError):
            (e_1 ^ e_2).get_coefs(3)

    def test_blade_coefs(self):
        """
        Various tests on several multivectors.
        """
        _g3d, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        m0 = 2 * e_1 + e_2 - e_3 + 3 * (e_1 ^ e_3) + (e_1 ^ e_3) + (e_2 ^ (3 * e_3))
        assert m0.blade_coefs([e_1]) == [2]
        assert m0.blade_coefs([e_2]) == [1]
        assert m0.blade_coefs([e_1, e_2]) == [2, 1]
        assert m0.blade_coefs([e_1 ^ e_3]) == [4]
        assert m0.blade_coefs([e_1 ^ e_3, e_2 ^ e_3]) == [4, 3]
        assert m0.blade_coefs([e_2 ^ e_3, e_1 ^ e_3]) == [3, 4]
        assert m0.blade_coefs([e_1, e_2 ^ e_3]) == [2, 3]

        a = sympy.Symbol('a')
        b = sympy.Symbol('b')
        m1 = a * e_1 + e_2 - e_3 + b * (e_1 ^ e_2)
        assert m1.blade_coefs([e_1]) == [a]
        assert m1.blade_coefs([e_2]) == [1]
        assert m1.blade_coefs([e_3]) == [-1]
        assert m1.blade_coefs([e_1, e_2, e_3]) == [a, 1, -1]
        assert m1.list() == [a, 1, -1]  # alias

        assert m1.blade_coefs([e_1 ^ e_2]) == [b]
        assert m1.blade_coefs([e_2 ^ e_3]) == [0]
        assert m1.blade_coefs([e_1 ^ e_3]) == [0]
        assert m1.blade_coefs([e_1 ^ e_2 ^ e_3]) == [0]

        # Invalid parameters
        pytest.raises(ValueError, lambda: m1.blade_coefs([e_1 + e_2]))
        pytest.raises(ValueError, lambda: m1.blade_coefs([e_2 ^ e_1]))
        pytest.raises(ValueError, lambda: m1.blade_coefs([e_1, e_2 ^ e_1]))
        pytest.raises(ValueError, lambda: m1.blade_coefs([a * e_1]))
        pytest.raises(ValueError, lambda: m1.blade_coefs([3 * e_3]))

    def test_rep_switching(self):
        # this ga has a non-diagonal metric
        _g3d, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        m0 = 2 * e_1 + e_2 - e_3 + 3 * (e_1 ^ e_3) + (e_1 ^ e_3) + (e_2 ^ (3 * e_3))
        m1 = (-4*(e_1 | e_3)-3*(e_2 | e_3))+2*e_1+e_2-e_3+4*e_1*e_3+3*e_2*e_3
        # m1 was chosen to make this true
        assert m0 == m1

        # all objects start off in blade rep
        assert m0.is_blade_rep

        # convert to base rep
        m0_base = m0.base_rep()
        assert m0.is_blade_rep  # original should not change
        assert not m0_base.is_blade_rep
        assert m0 == m0_base

        # convert back
        m0_base_blade = m0_base.blade_rep()
        assert not m0_base.is_blade_rep  # original should not change
        assert m0_base_blade.is_blade_rep
        assert m0 == m0_base_blade

    def test_construction(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        def check(x, expected_grades):
            assert x.grades == expected_grades
            assert x != 0

        # non-function symbol construction
        check(ga.mv('A', 'scalar'), [0])
        check(ga.mv('A', 'grade', 0), [0])
        check(ga.mv('A', 0), [0])
        check(ga.mv('A', 'vector'), [1])
        check(ga.mv('A', 'grade', 1), [1])
        check(ga.mv('A', 1), [1])
        check(ga.mv('A', 'bivector'), [2])
        check(ga.mv('A', 'grade2'), [2])
        check(ga.mv('A', 2), [2])
        check(ga.mv('A', 'pseudo'), [3])
        check(ga.mv('A', 'spinor'), [0, 2])
        check(ga.mv('A', 'even'), [0, 2])
        check(ga.mv('A', 'odd'), [1, 3])
        check(ga.mv('A', 'mv'), [0, 1, 2, 3])

        # value construction
        check(ga.mv([1, 2, 3], 'vector'), [1])

        # illegal arguments
        with pytest.raises(TypeError):
            ga.mv('A', 'vector', "too many arguments")
        with pytest.raises(TypeError):
            ga.mv('A', 'grade')  # too few arguments
        with pytest.raises(TypeError):
            ga.mv('A', 'grade', not_an_argument=True)  # invalid kwarg
        with pytest.raises(TypeError):
            ga.mv([1, 2, 3], 'vector', f=True)  # can't pass f with coefficients
        with pytest.raises(TypeError):
            ga.mv(e_1, 'even')  # Must be a string

    def test_abs(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3', g=[1, 1, 1])
        B = ga.mv('B', 'bivector')
        assert abs(B*B) == -(B*B).scalar()

    def test_hashable(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        d = {}
        d[e_1] = 1
        d[e_2] = 2
        assert d[e_1 + 0] == 1
        d[10] = 3  # note: not a multivector key!
        assert d[e_1 * 0 + 10] == 3

    def test_subs(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3', g=[1, 1, 1])
        B = ga.mv('B', 'bivector')
        B_inv = B.inv()

        B_mag = sympy.Symbol('|B|')

        # both of the sympy subs syntaxes work:
        assert (-B / B_mag**2).subs(B_mag, abs(B)) == B_inv
        assert (-B / B_mag**2).subs({B_mag: abs(B)}) == B_inv

    def test_sympify(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3', g=[1, 1, 1])

        # Letting this succeed silently and not return an Mv instance would be
        # dangerous.
        with pytest.raises(sympy.SympifyError):
            sympy.sympify(e_1)

        # this is fine
        sympy.sympify(e_1.obj)

    def test_arithmetic(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3', g=[1, 1, 1])
        one = ga.mv(sympy.S.One)

        # test that scalars are promoted to Mvs correctly
        assert e_1 + 1 == e_1 + one
        assert 1 + e_1 == one + e_1
        assert e_1 - 1 == e_1 - one
        assert 1 - e_1 == one - e_1

    @pytest.mark.parametrize('make_one', [
        lambda ga: 1,
        lambda ga: ga.mv(sympy.S.One),
        pytest.param(lambda ga: sympy.S.One, marks=pytest.mark.skipif(
            Version(sympy.__version__) < Version("1.6"),
            # until sympy/sympy@bec42df53cf2486d485065ddad1c31011a48bf3b
            reason="Cannot override < and > on sympy.Expr"
        ))
    ])
    def test_contraction(self, make_one):
        ga, e_1, e_2 = Ga.build('e*1|2', g=[1, 1])
        e12 = e_1 ^ e_2
        one = make_one(ga)

        assert (one < e12) == e12
        assert (e12 > one) == e12
        assert (e12 < one) == 0
        assert (one > e12) == 0

    def test_proj(self):
        g3coords = symbols('x y z', real=True)
        V = Ga('e', g=[1, 1, 1], coords=g3coords)

        u = V.mv("u", "vector")
        v = V.mv("v", "vector")
        w = V.mv("w", "vector")
        B = V.mv("B", "mv")

        assert proj(u, v) == v.project_in_blade(u)
        assert proj(w, v) + proj(w, u) == proj(w, u + v)
        assert proj(u, v) == (v | u) / u

        Vr = u ^ v
        assert proj(Vr, B) == ((B < Vr) * Vr.inv())
        assert proj(Vr, B) == ((B < Vr) < Vr.inv())

    def test_norm2(self):
        g3coords = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)
        A = g3.mv('A', 'mv')

        assert A.norm2() == A.rev().sp(A)
        assert A.norm2('+') == A.rev().sp(A)
        assert A.norm2('-') == A.rev().sp(A)

        assert norm2(A) == norm(A) * norm(A)
        assert norm2(A, '+') == norm(A) * norm(A)
        assert norm2(A, '-') == norm(A) * norm(A)

    def test_norm_nonneg(self):
        """Test that norm always returns a nonneg expression (issue 522)."""
        from sympy import Abs
        g3 = Ga('e', g=[1, 1, 1], coords=symbols('x y z', real=True))

        # scalar norm should include Abs and be nonneg
        s = g3.mv('s', 'scalar')
        s_norm = s.norm()
        assert isinstance(s_norm, Abs)
        assert s_norm.is_nonnegative

        # vector norm is already nonneg (sum of squares under sqrt)
        v = g3.mv('v', 'vector')
        assert v.norm().is_nonnegative

        # even-grade (bivector) norm is nonneg
        B = g3.mv('B', 'bivector')
        assert B.norm().is_nonnegative

        # non-Euclidean: symbolic vector where norm may not be provably nonneg
        g2mn = Ga('e', g=[1, -1], coords=symbols('x t', real=True))
        v = g2mn.mv('v', 'vector')
        assert v.norm(hint='+').is_nonnegative

    def test_mag2(self):
        g3coords = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)
        A = g3.mv('A', 'mv')

        assert mag2(A) == mag(A) * mag(A)

    def test_undual(self):
        # A not a multivector in undual(A).
        with pytest.raises(ValueError):
            undual(1)

    def test_g_invol(self):
        g3coords = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)
        A = g3.mv('A', 'mv')

        assert g_invol(A.even()) == A.even()
        assert g_invol(A.odd()) == -A.odd()

        with pytest.raises(ValueError):
            g_invol(1)

    def test_exp(self):
        g3coords = (x, y, z) = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)

        u = g3.mv("u", "vector")
        v = g3.mv("v", "vector")

        assert exp(u ^ v) == exp(-v ^ u)
        assert exp(x + y + z) == exp(z + y + x)

    def test_ccon(self):
        g3coords = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)

        A = g3.mv('A', 'mv')

        assert ccon(ccon(A)) == A
        assert ccon(A) == g_invol(rev(A))
        assert ccon(A) == rev(g_invol(A))

        # not a multivector in ccon(A)
        with pytest.raises(ValueError):
            ccon(1)

    def test_scalar(self):
        g3coords = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)

        A = g3.mv('A', 'mv')
        u = g3.mv("u", "vector")

        assert scalar(A) == A.scalar()
        assert scalar(u) == 0
        assert scalar(u * u) == qform(u)

        # not a multivector in scalar(A)
        with pytest.raises(ValueError):
            scalar(1)

    def test_sp(self):
        g3coords = symbols('x y z', real=True)
        g3 = Ga('e', g=[1, 1, 1], coords=g3coords)

        A = g3.mv('A', 'mv')
        B = g3.mv('B', 'mv')

        assert sp(A, B) == A.sp(B)
        assert sp(A, B) == (A*B).scalar()

        # not a multivector in sp(A, B)
        with pytest.raises(ValueError):
            sp(1, 1)

    @pytest.mark.parametrize('inv_func', [inv, shirokov_inverse, hitzer_inverse])
    @pytest.mark.parametrize('n', range(2, 5)) # tests are slow for larger n
    def test_inv(self, inv_func, n):
        gncoords = symbols(" ".join(f"x{i}" for i in range(n)), real=True)
        gn = Ga('e', g=[1] * n, coords=gncoords)

        # invert versors
        A = gn.mv('A', 'vector')
        Ainv = inv_func(A)
        assert A * Ainv == 1 + 0 * A

        # invert generic multivectors
        if inv_func in [shirokov_inverse, hitzer_inverse] and n<4:
            B = gn.mv('A', 'mv')
            Binv = inv_func(B)
            assert B * Binv == 1 + 0 * B

    def test_rtruediv(self):
        """Test that scalar/Mv works (issue 512)."""
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3')
        I = ga.I()

        # 1/I should be the inverse of I
        result = 1 / I
        assert result * I == ga.mv(1)

        # scalar / vector
        v = e_1 + e_2
        result = 2 / v
        assert result * v == ga.mv(2)

        # sympy scalar / Mv
        from sympy import S, symbols
        result = S(3) / I
        assert result * I == ga.mv(3)

        # symbolic numerator
        a = symbols('a')
        result = a / I
        assert result * I == ga.mv(a)

        # non-Euclidean metric (Minkowski-like)
        ga_m, e_t, e_x = Ga.build('e*t|x', g=[1, -1])
        v = e_t + e_x
        # v*v = 1 - 1 = 0 for null vector — skip (not invertible)
        v2 = e_t + 2 * e_x  # v2*v2 = 1 - 4 = -3 (invertible)
        result = 1 / v2
        assert result * v2 == ga_m.mv(1)

    # reproduce #537: is_blade() must handle null vectors and null blades

    def test_is_blade_zero_and_scalar(self):
        """Zero mv is not a blade; grade-0 scalars are."""
        ga, e1, e2 = Ga.build('e*1|2', g=[1, 1])
        assert ga.mv(S.Zero).is_blade() is False
        assert ga.mv(S.One).is_blade() is True

    def test_is_blade_vectors(self):
        """Grade-1: non-null and null vectors are both blades."""
        ga, e0, e1 = Ga.build('e*0|1', g=[1, -1])
        # non-null vectors
        assert e0.is_blade() is True
        assert e1.is_blade() is True
        # null vector b = e0+e1 (b*b = 0): is_versor() fails but is_blade() must pass
        b = e0 + e1
        assert b.is_zero() is False
        assert (b * b).is_zero()        # confirm null
        assert b.is_versor() is False   # no inverse
        assert b.is_blade() is True     # still a 1-blade

    def test_is_blade_bivectors(self):
        """Grade-2: non-null and null bivectors."""
        # non-null bivector via is_versor()
        ga, e1, e2, e3 = Ga.build('e*1|2|3', g=[1, 1, 1])
        assert (e1 ^ e2).is_blade() is True
        # null 2-blade: outer-product squaring test must detect it
        ga3, f0, f1, f2 = Ga.build('f*0|1|2', g=[1, -1, -1])
        B = (f0 + f1) ^ f2
        assert B.is_zero() is False
        assert B.is_versor() is False   # null constituent
        assert B.is_blade() is True     # but it IS a 2-blade

    def test_is_blade_trivector(self):
        """Grade-3: pseudoscalar of R^3 is a blade."""
        ga, e1, e2, e3 = Ga.build('e*1|2|3', g=[1, 1, 1])
        assert (e1 ^ e2 ^ e3).is_blade() is True

    @pytest.mark.slow
    def test_is_blade_grade3_known_limitation(self):
        """Grade >= 3: B^B=0 is not sufficient; is_blade() may give false positives.

        Dorst/Fontijne/Mann counterexample in R^6: the 3-vector
            X = e1^e2^e5 + e1^e3^e6 + e2^e4^e6 - e3^e4^e5
        satisfies X^X = 0 but is NOT a blade.  Our algorithm incorrectly
        returns True; see is_blade() docstring.
        """
        ga, e1, e2, e3, e4, e5, e6 = Ga.build('e*1|2|3|4|5|6', g=[1] * 6)
        X = (e1 ^ e2 ^ e5) + (e1 ^ e3 ^ e6) + (e2 ^ e4 ^ e6) - (e3 ^ e4 ^ e5)
        assert (X ^ X).is_zero()    # necessary condition satisfied → triggers false positive
        assert X.is_blade() is True  # known limitation: algorithm returns True

    def test_is_blade_non_simple_bivector(self):
        """Non-simple bivector: B^B != 0, so is_blade() correctly returns False."""
        ga, e1, e2, e3, e4 = Ga.build('e*1|2|3|4', g=[1, 1, 1, 1])
        B = (e1 ^ e2) + (e3 ^ e4)  # non-simple; B^B = 2*e1^e2^e3^e4 != 0
        assert not (B ^ B).is_zero()
        assert B.is_blade() is False

    def test_is_blade_non_homogeneous(self):
        """Non-grade-homogeneous mv is not a blade."""
        ga, e0, e1 = Ga.build('e*0|1', g=[1, -1])
        assert (e0 + (e0 ^ e1)).is_blade() is False

    def test_is_blade_result_cached(self):
        """is_blade() caches its result in blade_flg."""
        ga, e1, e2 = Ga.build('e*1|2', g=[1, 1])
        assert e1.is_blade() is True
        assert e1.blade_flg is True   # cache populated
        assert e1.is_blade() is True  # second call hits cached branch

    def test_reflect_in_null_blade_raises(self):
        """reflect_in_blade() must raise ValueError for null blades (#537)."""
        ga, f0, f1, f2 = Ga.build('f*0|1|2', g=[1, -1, -1])
        with pytest.raises(ValueError, match='null blade'):
            f0.reflect_in_blade((f0 + f1) ^ f2)

    def test_project_in_null_blade_raises(self):
        """project_in_blade() must raise ValueError for null blades (#537)."""
        ga, f0, f1, f2 = Ga.build('f*0|1|2', g=[1, -1, -1])
        with pytest.raises(ValueError, match='null blade'):
            f0.project_in_blade((f0 + f1) ^ f2)

