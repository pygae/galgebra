import unittest

from itertools import product
from sympy import simplify, Symbol, cos, sin, sqrt
from ga import Ga
from mv import Mv, cross


class TestChapter3(unittest.TestCase):

    def assertEquals(self, first, second, msg=None):
        """
        Compare two expressions are equals.
        """

        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        self.assertTrue(simplify(first - second) == 0, "%s == %s" % (first, second))


    def test6_1_4_1(self):
        """
        The geometric product for vectors on a basis.
        """
        R, e_1, e_2, e_3 = Ga.build('e*1|2|3', g='1 0 0, 0 1 0, 0 0 1')

        self.assertEquals(e_1 * e_1, 1)
        self.assertEquals(e_1 * e_2, e_1 ^ e_2)
        self.assertEquals(e_1 * e_3, e_1 ^ e_3)
        self.assertEquals(e_2 * e_1, -e_1 ^ e_2)
        self.assertEquals(e_2 * e_2, 1)
        self.assertEquals(e_2 * e_3, e_2 ^ e_3)
        self.assertEquals(e_3 * e_1, -e_1 ^ e_3)
        self.assertEquals(e_3 * e_2, -e_2 ^ e_3)
        self.assertEquals(e_3 * e_3, 1)

        e_12 = e_1 * e_2
        self.assertEquals(e_12 * e_12, -1)

        e_13 = e_1 * e_3
        self.assertEquals(e_13 * e_13, -1)

        e_23 = e_2 * e_3
        self.assertEquals(e_23 * e_23, -1)


    def test6_1_4_2(self):
        """
        The geometric product for vectors on a basis.
        """
        R, e_1, e_2 = Ga.build('e*1|2', g='1 0, 0 1')
        ONE = R.mv(1, 'scalar')
        e_12 = e_1 ^ e_2

        self.assertEquals(ONE * ONE, 1)
        self.assertEquals(ONE * e_1, e_1)
        self.assertEquals(ONE * e_2, e_2)
        self.assertEquals(ONE * e_12, e_12)

        self.assertEquals(e_1 * ONE, e_1)
        self.assertEquals(e_1 * e_1, 1)
        self.assertEquals(e_1 * e_2, e_12)
        self.assertEquals(e_1 * e_12, e_2)

        self.assertEquals(e_2 * ONE, e_2)
        self.assertEquals(e_2 * e_1, -e_12)
        self.assertEquals(e_2 * e_2, 1)
        self.assertEquals(e_2 * e_12, -e_1)

        self.assertEquals(e_12 * ONE, e_12)
        self.assertEquals(e_12 * e_1, -e_2)
        self.assertEquals(e_12 * e_2, e_1)
        self.assertEquals(e_12 * e_12, -1)

        a_1 = Symbol('a_1')
        a_2 = Symbol('a_2')
        a = R.mv((a_1, a_2), 'vector')

        b_1 = Symbol('b_1')
        b_2 = Symbol('b_2')
        b = R.mv((b_1, b_2), 'vector')

        self.assertEquals(a * b, a_1 * b_1 + a_2 * b_2 + (a_1 * b_2 - a_2 * b_1) * (e_1 ^ e_2))


    def test6_3_1(self):
        """
        The subspace products from symmetry.
        """
        R = Ga('e*1|2|3')

        a = R.mv('a', 'vector')
        B_blades = [R.mv('B%d' % k, k, 'grade') for k in range(R.n + 1)]

        def hat(M):
            M_grades = R.grade_decomposition(M).keys()
            self.assertEquals(len(M_grades), 1)
            return ((-1) ** M_grades[0]) * M

        for B in B_blades:
            self.assertEquals(a ^ B, 0.5 * (a * B + hat(B) * a))
            self.assertEquals(B ^ a, 0.5 * (B * a + a * hat(B)))
            self.assertEquals(a < B, 0.5 * (a * B - hat(B) * a))
            self.assertEquals(B > a, 0.5 * (B * a - a * hat(B)))
            self.assertEquals(a * B, (a < B) + (a ^ B))

        for B in B_blades:
            self.assertEquals(hat(B) * a, (hat(B) > a) + (hat(B) ^ a))
            self.assertEquals(hat(B) * a, -(a < B) + (a ^ B))
            self.assertEquals(hat(B) * a, a * B - 2 * (a < B))
            self.assertEquals(hat(B) * a, -(a * B) + 2 * (a ^ B))

        for B in B_blades:
            self.assertEquals(a * B, hat(B) * a + 2 * (a < B))
            self.assertEquals(a * B, -hat(B) * a + 2 * (a ^ B))

        M = R.mv('m', 'mv')
        self.assertEquals(a * M, (a < M) + (a ^ M))
