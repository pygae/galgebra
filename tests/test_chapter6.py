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


    def test6_3_2(self):
        """
        The subspace products as selected grades.
        """

        R = Ga('e*1|2|3')

        A_blades = [R.mv('A%d' % k, k, 'grade') for k in range(R.n + 1)]
        B_blades = [R.mv('B%d' % l, l, 'grade') for l in range(R.n + 1)]

        def grade(M):
            M_grades = R.grade_decomposition(M).keys()
            self.assertEquals(len(M_grades), 1)
            return M_grades[0]

        for A, B in product(A_blades, B_blades):
            k = grade(A)
            l = grade(B)
            self.assertEquals(A ^ B, (A * B).get_grade(k + l))
            self.assertEquals(A < B, 0 if k > l else (A * B).get_grade(l - k))
            self.assertEquals(A > B, 0 if l > k else (A * B).get_grade(k - l))
            # TODO: scalar product


    def test6_6_2_3(self):
        """
        The outer product can be defined as the completely antisymmetric summed average of all permutations
        of the geometric product of its factors, with a sign for each term depending on oddness or evenness of
        the permutation. Derive x ^ y ^ z = 1/3! * (xyz - yxz + yzx - zyx + zxy - xzy)
        """

        R = Ga('e*1|2|3')

        x = R.mv('x', 'vector')
        y = R.mv('y', 'vector')
        z = R.mv('z', 'vector')

        def hat(M):
            M_grades = R.grade_decomposition(M).keys()
            self.assertEquals(len(M_grades), 1)
            return ((-1) ** M_grades[0]) * M
        
        self.assertEquals(x * y * z, ((x < y) + (x ^ y)) * z)
        self.assertEquals(x * y * z, ((x < y) * z) + ((x ^ y) * z))
        self.assertEquals(x * y * z, ((x < y) * z) - (z < hat(x ^ y)) + (z ^ hat(x ^ y)))
        self.assertEquals(hat(x ^ y), x ^ y)
        self.assertEquals(x * y * z, ((x < y) * z) - (z < (x ^ y)) + (z ^ x ^ y))
        self.assertEquals(x * y * z, ((x < y) * z) - (z < (x ^ y)) + (x ^ y ^ z))

        self.assertEquals(y * x * z, ((y < x) + (y ^ x)) * z)
        self.assertEquals(y * x * z, ((y < x) * z) + ((y ^ x) * z))
        self.assertEquals(y * x * z, ((y < x) * z) - (z < hat(y ^ x)) + (z ^ hat(y ^ x)))
        self.assertEquals(hat(y ^ x), y ^ x)
        self.assertEquals(y * x * z, ((y < x) * z) - (z < (y ^ x)) + (z ^ y ^ x))
        self.assertEquals(y * x * z, ((x < y) * z) + (z < (x ^ y)) - (x ^ y ^ z))

        self.assertEquals(y * z * x, ((y < z) + (y ^ z)) * x)
        self.assertEquals(y * z * x, ((y < z) * x) - (x < hat(y ^ z)) + (x ^ hat(y ^ z)))
        self.assertEquals(y * z * x, ((y < z) * x) - (x < (y ^ z)) + (x ^ y ^ z))

        self.assertEquals(z * y * x, ((z < y) + (z ^ y)) * x)
        self.assertEquals(z * y * x, ((z < y) * x) - (x < hat(z ^ y)) + (x ^ hat(z ^ y)))
        self.assertEquals(z * y * x, ((z < y) * x) - (x < (z ^ y)) - (x ^ y ^ z))
        self.assertEquals(z * y * x, ((y < z) * x) + (x < (y ^ z)) - (x ^ y ^ z))

        self.assertEquals(z * x * y, ((z < x) + (z ^ x)) * y)
        self.assertEquals(z * x * y, ((z < x) * y) + ((z ^ x) * y))
        self.assertEquals(z * x * y, ((z < x) * y) - (y < hat(z ^ x)) + (y ^ hat(z ^ x)))
        self.assertEquals(z * x * y, ((z < x) * y) - (y < (z ^ x)) + (y ^ z ^ x))
        self.assertEquals(z * x * y, ((z < x) * y) - (y < (z ^ x)) + (x ^ y ^ z))

        self.assertEquals(x * z * y, ((x < z) + (x ^ z)) * y)
        self.assertEquals(x * z * y, ((x < z) * y) + ((x ^ z) * y))
        self.assertEquals(x * z * y, ((x < z) * y) - (y < hat(x ^ z)) + (y ^ hat(x ^ z)))
        self.assertEquals(x * z * y, ((x < z) * y) - (y < (x ^ z)) + (y ^ x ^ z))
        self.assertEquals(x * z * y, ((z < x) * y) + (y < (z ^ x)) - (x ^ y ^ z))

        self.assertEquals(x * y * z - y * x * z, 2 * ((x ^ y ^ z) - (z < (x ^ y))))
        self.assertEquals(y * z * x - z * y * x, 2 * ((x ^ y ^ z) - (x < (y ^ z))))
        self.assertEquals(z * x * y - x * z * y, 2 * ((x ^ y ^ z) - (y < (z ^ x))))

        self.assertEquals(z < (x ^ y), ((z < x) ^ y) - (x ^ (z < y)))
        self.assertEquals(x < (y ^ z), ((x < y) ^ z) - (y ^ (x < z)))
        self.assertEquals(y < (z ^ x), ((y < z) ^ x) - (z ^ (y < x)))

        self.assertEquals(((z < x) ^ y) - (x ^ (z < y)) + ((x < y) ^ z) - (y ^ (x < z)) + ((y < z) ^ x) - (z ^ (y < x)), 0)

        self.assertEquals(x * y * z - y * x * z + y * z * x - z * y * x + z * x * y - x * z * y, 6 * (x ^ y ^ z))

