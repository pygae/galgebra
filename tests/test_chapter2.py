import sys
sys.path.append('../galgebra')

import unittest

from sympy import Symbol, solve
from ga import Ga

class TestChapter2(unittest.TestCase):

    def test_2_12_1_1(self):
        """
        Compute the outer products of the following 3-space expressions,
        giving the result relative to the basis {1, e_1, e_2, e_3, e_1^e_2, e_1^e_3, e_2^e_3, e_1^e_2^e_3}.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')
        self.assertTrue((e_1 + e_2) ^ (e_1 + e_3) == (-e_1 ^ e_2) + (e_1 ^ e_3) + (e_2 ^ e_3))
        self.assertTrue((e_1 + e_2 + e_3) ^ (2*e_1) == -2*(e_1 ^ e_2) - 2*(e_1 ^ e_3))
        self.assertTrue((e_1 - e_2) ^ (e_1 - e_3) == (e_1 ^ e_2) - (e_1 ^ e_3) + (e_2 ^ e_3))
        self.assertTrue((e_1 + e_2) ^ (0.5*e_1 + 2*e_2 + 3*e_3) == 1.5*(e_1 ^ e_2) + 3*(e_1 ^ e_3) + 3*(e_2 ^ e_3))
        self.assertTrue((e_1 ^ e_2) ^ (e_1 + e_3) == (e_1 ^ e_2 ^ e_3))
        self.assertTrue((e_1 + e_2) ^ ((e_1 ^ e_2) + (e_2 ^ e_3)) == (e_1 ^ e_2 ^ e_3))


    def test2_12_1_2(self):
        """
        Given the 2-blade B = e_1 ^ (e_2 - e_3) that represents a plane,
        determine if each of the following vectors lies in that plane.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')
        B = e_1 ^ (e_2 - e_3)
        self.assertTrue(e_1 ^ B == 0)
        self.assertFalse((e_1 + e_2) ^ B == 0)
        self.assertFalse((e_1 + e_2 + e_3) ^ B == 0)
        self.assertTrue((2*e_1 - e_2 + e_3) ^ B == 0)


    def test2_12_1_3(self):
        """
        What is the area of the parallelogram spanned by the vectors a = e_1 + 2*e_2
        and b = -e_1 - e_2 (relative to the area of e_1 ^ e_2) ?
        """
        (_g3d, e_1, e_2, _e_3) = Ga.build('e*1|2|3')
        a = e_1 + 2*e_2
        b = -e_1 - e_2
        B = a ^ b
        self.assertTrue(B == 1 * (e_1 ^ e_2))


    def test2_12_1_4(self):
        """
        Compute the intersection of the non-homogeneous line L with position vector e_1
        and direction vector e_2, and the line M with position vector e_2 and direction vector
        (e_1 + e_2), using 2-blades.
        """
        (_g2d, e_1, e_2) = Ga.build('e*1|2')

        # x is defined in the basis {e_1, e_2}
        a = Symbol('a')
        b = Symbol('b')
        x = a * e_1 + b * e_2

        # x lies on L and M (eq. L == 0 and M == 0)
        L = (x ^ e_2) - (e_1 ^ e_2)
        M = (x ^ (e_1 + e_2)) - (e_2 ^ (e_1 + e_2))

        # Solve the linear system
        R = solve([L, M], a, b)

        # Replace symbols
        x = x.subs(R)

        self.assertTrue(x == e_1 + 2*e_2)


if __name__ == '__main__':

    unittest.main()

