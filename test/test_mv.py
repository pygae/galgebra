from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import unittest
from sympy import Symbol
from galgebra.ga import Ga

class TestMv(unittest.TestCase):

    def test_is_base(self):
        """
        Various tests on several multivectors.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')

        self.assertTrue((e_1).is_base())
        self.assertTrue((e_2).is_base())
        self.assertTrue((e_3).is_base())
        self.assertTrue((e_1 ^ e_2).is_base())
        self.assertTrue((e_2 ^ e_3).is_base())
        self.assertTrue((e_1 ^ e_3).is_base())
        self.assertTrue((e_1 ^ e_2 ^ e_3).is_base())

        self.assertFalse((2*e_1).is_base())
        self.assertFalse((e_1 + e_2).is_base())
        self.assertFalse((e_3 * 4).is_base())
        self.assertFalse(((3 * e_1) ^ e_2).is_base())
        self.assertFalse((2 * (e_2 ^ e_3)).is_base())
        self.assertFalse((e_3 ^ e_1).is_base())
        self.assertFalse((e_2 ^ e_1 ^ e_3).is_base())


    def test_blade_coefs(self):
        """
        Various tests on several multivectors.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')

        m0 =  2 * e_1 + e_2 - e_3 + 3 * (e_1 ^ e_3) + (e_1 ^ e_3) + (e_2 ^ (3 * e_3))
        self.assertTrue(m0.blade_coefs([e_1]) == [2])
        self.assertTrue(m0.blade_coefs([e_2]) == [1])
        self.assertTrue(m0.blade_coefs([e_1, e_2]) == [2, 1])
        self.assertTrue(m0.blade_coefs([e_1 ^ e_3]) == [4])
        self.assertTrue(m0.blade_coefs([e_1 ^ e_3, e_2 ^ e_3]) == [4, 3])
        self.assertTrue(m0.blade_coefs([e_2 ^ e_3, e_1 ^ e_3]) == [3, 4])
        self.assertTrue(m0.blade_coefs([e_1, e_2 ^ e_3]) == [2, 3])

        a = Symbol('a')
        b = Symbol('b')
        m1 = a * e_1 + e_2 - e_3 + b * (e_1 ^ e_2)
        self.assertTrue(m1.blade_coefs([e_1]) == [a])
        self.assertTrue(m1.blade_coefs([e_2]) == [1])
        self.assertTrue(m1.blade_coefs([e_3]) == [-1])
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
