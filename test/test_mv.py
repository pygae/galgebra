import unittest

import pytest
import sympy
from galgebra.ga import Ga

class TestMv(unittest.TestCase):

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

        self.assertTrue((e_1).is_base())
        self.assertTrue((e_2).is_base())
        self.assertTrue((e_3).is_base())
        self.assertTrue((e_1 ^ e_2).is_base())
        self.assertTrue((e_2 ^ e_3).is_base())
        self.assertTrue((e_1 ^ e_3).is_base())
        self.assertTrue((e_1 ^ e_2 ^ e_3).is_base())

        self.assertFalse((2 * e_1).is_base())
        self.assertFalse((e_1 + e_2).is_base())
        self.assertFalse((e_3 * 4).is_base())
        self.assertFalse(((3 * e_1) ^ e_2).is_base())
        self.assertFalse((2 * (e_2 ^ e_3)).is_base())
        self.assertFalse((e_3 ^ e_1).is_base())
        self.assertFalse((e_2 ^ e_1 ^ e_3).is_base())

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

        m0 =  2 * e_1 + e_2 - e_3 + 3 * (e_1 ^ e_3) + (e_1 ^ e_3) + (e_2 ^ (3 * e_3))
        self.assertTrue(m0.blade_coefs([e_1]) == [2])
        self.assertTrue(m0.blade_coefs([e_2]) == [1])
        self.assertTrue(m0.blade_coefs([e_1, e_2]) == [2, 1])
        self.assertTrue(m0.blade_coefs([e_1 ^ e_3]) == [4])
        self.assertTrue(m0.blade_coefs([e_1 ^ e_3, e_2 ^ e_3]) == [4, 3])
        self.assertTrue(m0.blade_coefs([e_2 ^ e_3, e_1 ^ e_3]) == [3, 4])
        self.assertTrue(m0.blade_coefs([e_1, e_2 ^ e_3]) == [2, 3])

        a = sympy.Symbol('a')
        b = sympy.Symbol('b')
        m1 = a * e_1 + e_2 - e_3 + b * (e_1 ^ e_2)
        self.assertTrue(m1.blade_coefs([e_1]) == [a])
        self.assertTrue(m1.blade_coefs([e_2]) == [1])
        self.assertTrue(m1.blade_coefs([e_3]) == [-1])
        self.assertTrue(m1.blade_coefs([e_1, e_2, e_3]) == [a, 1, -1])
        self.assertTrue(m1.list() == [a, 1, -1])  # alias

        self.assertTrue(m1.blade_coefs([e_1 ^ e_2]) == [b])
        self.assertTrue(m1.blade_coefs([e_2 ^ e_3]) == [0])
        self.assertTrue(m1.blade_coefs([e_1 ^ e_3]) == [0])
        self.assertTrue(m1.blade_coefs([e_1 ^ e_2 ^ e_3]) == [0])

        # Invalid parameters
        self.assertRaises(ValueError, lambda: m1.blade_coefs([e_1 + e_2]))
        self.assertRaises(ValueError, lambda: m1.blade_coefs([e_2 ^ e_1]))
        self.assertRaises(ValueError, lambda: m1.blade_coefs([e_1, e_2 ^ e_1]))
        self.assertRaises(ValueError, lambda: m1.blade_coefs([a * e_1]))
        self.assertRaises(ValueError, lambda: m1.blade_coefs([3 * e_3]))

    def test_rep_switching(self):
        # this ga has a non-diagonal metric
        _g3d, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        m0 =  2 * e_1 + e_2 - e_3 + 3 * (e_1 ^ e_3) + (e_1 ^ e_3) + (e_2 ^ (3 * e_3))
        m1 = (-4*(e_1 | e_3)-3*(e_2 | e_3))+2*e_1+e_2-e_3+4*e_1*e_3+3*e_2*e_3
        # m1 was chosen to make this true
        self.assertEqual(m0, m1)

        # all objects start off in blade rep
        self.assertTrue(m0.is_blade_rep)

        # convert to base rep
        m0_base = m0.base_rep()
        self.assertTrue(m0.is_blade_rep)  # original should not change
        self.assertFalse(m0_base.is_blade_rep)
        self.assertEqual(m0, m0_base)

        # convert back
        m0_base_blade = m0_base.blade_rep()
        self.assertFalse(m0_base.is_blade_rep)  # original should not change
        self.assertTrue(m0_base_blade.is_blade_rep)
        self.assertEqual(m0, m0_base_blade)

    def test_construction(self):
        ga, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        def check(x, expected_grades):
            self.assertEqual(x.grades, expected_grades)
            self.assertNotEqual(x, 0)

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
        with self.assertRaises(TypeError):
            ga.mv('A', 'vector', "too many arguments")
        with self.assertRaises(TypeError):
            ga.mv('A', 'grade')  # too few arguments
        with self.assertRaises(TypeError):
            ga.mv('A', 'grade', not_an_argument=True)  # invalid kwarg
        with self.assertRaises(TypeError):
            ga.mv([1, 2, 3], 'vector', f=True)  # can't pass f with coefficients

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
