from .test_utils import TestCase

from itertools import product
from sympy import Symbol, solve_poly_system
from galgebra.ga import Ga


class TestChapter6(TestCase):

    def test6_1_4_1(self):
        """
        The geometric product for vectors on a basis.
        """
        GA, e_1, e_2, e_3 = Ga.build('e*1|2|3', g='1 0 0, 0 1 0, 0 0 1')

        self.assertEqual(e_1 * e_1, 1)
        self.assertEqual(e_1 * e_2, e_1 ^ e_2)
        self.assertEqual(e_1 * e_3, e_1 ^ e_3)
        self.assertEqual(e_2 * e_1, -e_1 ^ e_2)
        self.assertEqual(e_2 * e_2, 1)
        self.assertEqual(e_2 * e_3, e_2 ^ e_3)
        self.assertEqual(e_3 * e_1, -e_1 ^ e_3)
        self.assertEqual(e_3 * e_2, -e_2 ^ e_3)
        self.assertEqual(e_3 * e_3, 1)

        e_12 = e_1 * e_2
        self.assertEqual(e_12 * e_12, -1)

        e_13 = e_1 * e_3
        self.assertEqual(e_13 * e_13, -1)

        e_23 = e_2 * e_3
        self.assertEqual(e_23 * e_23, -1)

    def test6_1_4_2(self):
        """
        The geometric product for vectors on a basis.
        """
        GA, e_1, e_2 = Ga.build('e*1|2', g='1 0, 0 1')
        ONE = GA.mv(1, 'scalar')
        e_12 = e_1 ^ e_2

        self.assertEqual(ONE * ONE, 1)
        self.assertEqual(ONE * e_1, e_1)
        self.assertEqual(ONE * e_2, e_2)
        self.assertEqual(ONE * e_12, e_12)

        self.assertEqual(e_1 * ONE, e_1)
        self.assertEqual(e_1 * e_1, 1)
        self.assertEqual(e_1 * e_2, e_12)
        self.assertEqual(e_1 * e_12, e_2)

        self.assertEqual(e_2 * ONE, e_2)
        self.assertEqual(e_2 * e_1, -e_12)
        self.assertEqual(e_2 * e_2, 1)
        self.assertEqual(e_2 * e_12, -e_1)

        self.assertEqual(e_12 * ONE, e_12)
        self.assertEqual(e_12 * e_1, -e_2)
        self.assertEqual(e_12 * e_2, e_1)
        self.assertEqual(e_12 * e_12, -1)

        a_1 = Symbol('a_1', real=True)
        a_2 = Symbol('a_2', real=True)
        a = GA.mv((a_1, a_2), 'vector')

        b_1 = Symbol('b_1', real=True)
        b_2 = Symbol('b_2', real=True)
        b = GA.mv((b_1, b_2), 'vector')

        self.assertEqual(a * b, a_1 * b_1 + a_2 * b_2 + (a_1 * b_2 - a_2 * b_1) * (e_1 ^ e_2))

    def test6_3_1(self):
        """
        The subspace products from symmetry.
        """
        GA = Ga('e*1|2|3')
        a = GA.mv('a', 'vector')
        B_grade_and_blades = [(k, GA.mv('B%d' % k, 'blade', k)) for k in range(GA.n + 1)]

        def hat(k, M):
            return ((-1) ** k) * M

        for k, B in B_grade_and_blades:
            self.assertEqual(a ^ B, 0.5 * (a * B + hat(k, B) * a))
            self.assertEqual(B ^ a, 0.5 * (B * a + a * hat(k, B)))
            self.assertEqual(a < B, 0.5 * (a * B - hat(k, B) * a))
            self.assertEqual(B > a, 0.5 * (B * a - a * hat(k, B)))
            self.assertEqual(a * B, (a < B) + (a ^ B))

        for k, B in B_grade_and_blades:
            self.assertEqual(hat(k, B) * a, (hat(k, B) > a) + (hat(k, B) ^ a))
            self.assertEqual(hat(k, B) * a, -(a < B) + (a ^ B))
            self.assertEqual(hat(k, B) * a, a * B - 2 * (a < B))
            self.assertEqual(hat(k, B) * a, -(a * B) + 2 * (a ^ B))

        for k, B in B_grade_and_blades:
            self.assertEqual(a * B, hat(k, B) * a + 2 * (a < B))
            self.assertEqual(a * B, -hat(k, B) * a + 2 * (a ^ B))

        M = GA.mv('m', 'mv')
        self.assertEqual(a * M, (a < M) + (a ^ M))

    def test6_3_2(self):
        """
        The subspace products as selected grades.
        """
        GA = Ga('e*1|2|3')

        # This test should work for blades and other graded multivectors
        A_grade_and_blades = [(k, GA.mv('A%d' % k, 'grade', k)) for k in range(GA.n + 1)]
        B_grade_and_blades = [(l, GA.mv('B%d' % l, 'grade', l)) for l in range(GA.n + 1)]

        for (k, A), (l, B) in product(A_grade_and_blades, B_grade_and_blades):
            self.assertEqual(A ^ B, (A * B).get_grade(k + l))
            self.assertEqual(A < B, 0 if k > l else (A * B).get_grade(l - k))
            self.assertEqual(A > B, 0 if l > k else (A * B).get_grade(k - l))
            # TODO: scalar product

    def test6_6_2_3(self):
        """
        The outer product can be defined as the completely antisymmetric summed average of all permutations
        of the geometric product of its factors, with a sign for each term depending on oddness or evenness of
        the permutation. Derive x ^ y ^ z = 1/3! * (xyz - yxz + yzx - zyx + zxy - xzy)
        """
        GA = Ga('e*1|2|3')
        x = GA.mv('x', 'vector')
        y = GA.mv('y', 'vector')
        z = GA.mv('z', 'vector')

        def hat(M):
            M_grades = list(GA.grade_decomposition(M).keys())
            self.assertEqual(len(M_grades), 1)
            return ((-1) ** M_grades[0]) * M

        self.assertEqual(x * y * z, ((x < y) + (x ^ y)) * z)
        self.assertEqual(x * y * z, ((x < y) * z) + ((x ^ y) * z))
        self.assertEqual(x * y * z, ((x < y) * z) - (z < hat(x ^ y)) + (z ^ hat(x ^ y)))
        self.assertEqual(hat(x ^ y), x ^ y)
        self.assertEqual(x * y * z, ((x < y) * z) - (z < (x ^ y)) + (z ^ x ^ y))
        self.assertEqual(x * y * z, ((x < y) * z) - (z < (x ^ y)) + (x ^ y ^ z))

        self.assertEqual(y * x * z, ((y < x) + (y ^ x)) * z)
        self.assertEqual(y * x * z, ((y < x) * z) + ((y ^ x) * z))
        self.assertEqual(y * x * z, ((y < x) * z) - (z < hat(y ^ x)) + (z ^ hat(y ^ x)))
        self.assertEqual(hat(y ^ x), y ^ x)
        self.assertEqual(y * x * z, ((y < x) * z) - (z < (y ^ x)) + (z ^ y ^ x))
        self.assertEqual(y * x * z, ((x < y) * z) + (z < (x ^ y)) - (x ^ y ^ z))

        self.assertEqual(y * z * x, ((y < z) + (y ^ z)) * x)
        self.assertEqual(y * z * x, ((y < z) * x) - (x < hat(y ^ z)) + (x ^ hat(y ^ z)))
        self.assertEqual(y * z * x, ((y < z) * x) - (x < (y ^ z)) + (x ^ y ^ z))

        self.assertEqual(z * y * x, ((z < y) + (z ^ y)) * x)
        self.assertEqual(z * y * x, ((z < y) * x) - (x < hat(z ^ y)) + (x ^ hat(z ^ y)))
        self.assertEqual(z * y * x, ((z < y) * x) - (x < (z ^ y)) - (x ^ y ^ z))
        self.assertEqual(z * y * x, ((y < z) * x) + (x < (y ^ z)) - (x ^ y ^ z))

        self.assertEqual(z * x * y, ((z < x) + (z ^ x)) * y)
        self.assertEqual(z * x * y, ((z < x) * y) + ((z ^ x) * y))
        self.assertEqual(z * x * y, ((z < x) * y) - (y < hat(z ^ x)) + (y ^ hat(z ^ x)))
        self.assertEqual(z * x * y, ((z < x) * y) - (y < (z ^ x)) + (y ^ z ^ x))
        self.assertEqual(z * x * y, ((z < x) * y) - (y < (z ^ x)) + (x ^ y ^ z))

        self.assertEqual(x * z * y, ((x < z) + (x ^ z)) * y)
        self.assertEqual(x * z * y, ((x < z) * y) + ((x ^ z) * y))
        self.assertEqual(x * z * y, ((x < z) * y) - (y < hat(x ^ z)) + (y ^ hat(x ^ z)))
        self.assertEqual(x * z * y, ((x < z) * y) - (y < (x ^ z)) + (y ^ x ^ z))
        self.assertEqual(x * z * y, ((z < x) * y) + (y < (z ^ x)) - (x ^ y ^ z))

        self.assertEqual(x * y * z - y * x * z, 2 * ((x ^ y ^ z) - (z < (x ^ y))))
        self.assertEqual(y * z * x - z * y * x, 2 * ((x ^ y ^ z) - (x < (y ^ z))))
        self.assertEqual(z * x * y - x * z * y, 2 * ((x ^ y ^ z) - (y < (z ^ x))))

        self.assertEqual(z < (x ^ y), ((z < x) ^ y) - (x ^ (z < y)))
        self.assertEqual(x < (y ^ z), ((x < y) ^ z) - (y ^ (x < z)))
        self.assertEqual(y < (z ^ x), ((y < z) ^ x) - (z ^ (y < x)))

        self.assertEqual(((z < x) ^ y) - (x ^ (z < y)) + ((x < y) ^ z) - (y ^ (x < z)) + ((y < z) ^ x) - (z ^ (y < x)),
                         0)

        self.assertEqual(x * y * z - y * x * z + y * z * x - z * y * x + z * x * y - x * z * y, 6 * (x ^ y ^ z))

    def test6_6_2_4(self):
        """
        The parts of a certain grade of a geometric product of blades are not necessarily blades.
        Show that in a 4D space with orthogonal basis, a counterexample is the grade 2 of the product
        e1 * (e1 + e2) * (e2 + e3) * (e1 + e4).
        """
        GA, e_1, e_2, e_3, e_4 = Ga.build('e*1|2|3|4', g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')
        R = (e_1 * (e_1 + e_2) * (e_2 + e_3) * (e_1 + e_4)).get_grade(2)

        # A 2-blade
        a_1 = Symbol('a_1', real=True)
        a_2 = Symbol('a_2', real=True)
        a_3 = Symbol('a_3', real=True)
        a_4 = Symbol('a_4', real=True)
        a = GA.mv((a_1, a_2, a_3, a_4), 'vector')

        b_1 = Symbol('b_1', real=True)
        b_2 = Symbol('b_2', real=True)
        b_3 = Symbol('b_3', real=True)
        b_4 = Symbol('b_4', real=True)
        b = GA.mv((b_1, b_2, b_3, b_4), 'vector')
        S = a ^ b

        # Try to solve the system and show there is no solution
        system = [
            S_coef - R_coef for S_coef, R_coef in zip(S.blade_coefs(), R.blade_coefs())
        ]

        unknowns = [
            a_1, a_2, a_3, a_4, b_1, b_2, b_3, b_4
        ]

        # TODO: use solve if sympy fix it
        result = solve_poly_system(system, unknowns)
        self.assertTrue(result is None)

    def test6_6_2_6(self):
        """
        Prove (X ^ A) * B = X * (A < B) using the grade based definition of ^, * and <.
        """
        GA = Ga('e*1|2|3')
        X_grade_blades = [(j, GA.mv('X', 'blade', j)) for j in range(GA.n + 1)]
        A_grade_blades = [(k, GA.mv('A', 'blade', k)) for k in range(GA.n + 1)]
        B_grade_blades = [(l, GA.mv('B', 'blade', l)) for l in range(GA.n + 1)]

        for (j, X), (k, A), (l, B) in product(X_grade_blades, A_grade_blades, B_grade_blades):
            self.assertEqual((X * A).get_grade(j + k), X ^ A)
            self.assertEqual((A * B).get_grade(l - k), 0 if k > l else A < B)
            self.assertTrue(((X * A).get_grade(j + k) * B).get_grade(0) != 0 if j + k == l else True)
            self.assertTrue((X * (A * B).get_grade(l - k)).get_grade(0) != 0 if j == l - k else True)
            self.assertEqual(((X * A).get_grade(j + k) * B).get_grade(0), (X * (A * B).get_grade(l - k)).get_grade(0))

    def test6_6_2_7(self):
        """
        In the formula (x < (1/A)) * A, show we can replace the geometric product by a contraction, so that
        it is in fact the projection (x < (1/A)) < A.
        """

        GA_list = [
            Ga('e*1|2', g='1 0, 0 1'),
            Ga('e*1|2|3', g='1 0 0, 0 1 0, 0 0 1'),
            Ga('e*1|2|3|4', g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1'),
            # Ga('e*1|2|3|4|5', g='1 0 0 0 0, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, 0 0 0 0 1'),
        ]

        for GA in GA_list:
            x = GA.mv('x', 'vector')
            A_grade_and_blades = [(k, GA.mv('A', 'blade', k)) for k in range(GA.n + 1)]
            for k, A in A_grade_and_blades:
                self.assertEqual((x < A.inv()) * A, (x < A.inv()) < A)

    def test6_6_2_9(self):
        """
        In a 4D space with orthonormal basis {e1, e2, e3, e4}, project the 2-blade X = (e1 + e2) ^ (e3 + e4) onto
        the 2-blade A = (e1 ^ e3). Then determine the rejection as the difference of X and its projection. Show
        that is not a blade.
        """
        GA, e_1, e_2, e_3, e_4 = Ga.build('e*1|2|3|4', g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')
        X = (e_1 + e_2) ^ (e_3 + e_4)
        A = (e_1 ^ e_3)

        # projection of X onto B
        def proj(X, B):
            return (X < B.inv()) < B

        P = proj(X, A)
        R = X - P

        # A 2-blade
        a_1 = Symbol('a_1', real=True)
        a_2 = Symbol('a_2', real=True)
        a_3 = Symbol('a_3', real=True)
        a_4 = Symbol('a_4', real=True)
        a = GA.mv((a_1, a_2, a_3, a_4), 'vector')

        b_1 = Symbol('b_1', real=True)
        b_2 = Symbol('b_2', real=True)
        b_3 = Symbol('b_3', real=True)
        b_4 = Symbol('b_4', real=True)
        b = GA.mv((b_1, b_2, b_3, b_4), 'vector')
        S = a ^ b

        # Try to solve the system and show there is no solution
        system = [
            S_coef - R_coef for S_coef, R_coef in zip(S.blade_coefs(), R.blade_coefs())
        ]

        unknowns = [
            a_1, a_2, a_3, a_4, b_1, b_2, b_3, b_4
        ]

        # TODO: use solve if sympy fix it
        result = solve_poly_system(system, unknowns)
        self.assertTrue(result is None)
