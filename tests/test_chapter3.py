import unittest

from itertools import product
from sympy import simplify, Symbol
from ga import Ga
from mv import Mv


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


    def test_3_2_2(self):
        """
        Computing the contraction explicitly.
        """

        Ga.dual_mode("Iinv+")

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]
        B_blades = [R.mv('B', i, 'grade') for i in range(R.n + 1)]
        C_blades = [R.mv('C', i, 'grade') for i in range(R.n + 1)]

        # scalar and blades of various grades
        A = A_blades[0]
        for B in B_blades:
            self.assertEquals(A < B, A * B)

        A = A_blades[0]
        for B in B_blades:
            self.assertEquals(B < A, 0 if B.pure_grade() > 0 else A * B)

        # vectors
        A = A_blades[1]
        B = B_blades[1]
        self.assertEquals(A < B, A | B)

        # vector and the outer product of 2 blades of various grades (scalars, vectors, 2-vectors...)
        A = A_blades[1]
        for B, C in product(B_blades, C_blades):
            self.assertEquals(A < (B ^ C), ((A < B) ^ C) + (-1)**B.pure_grade() * (B ^ (A < C)))

        # vector and the outer product of 2 blades of various grades (scalars, vectors, 2-vectors...)
        for A, B, C in product(A_blades, B_blades, C_blades):
            self.assertEquals((A ^ B) < C, A < (B < C))

        # distributive properties
        for A, B, C in product(A_blades, B_blades, C_blades):
            self.assertEquals((A + B) < C, (A < C) + (B < C))

        for A, B, C in product(A_blades, B_blades, C_blades):
            self.assertEquals(A < (B + C), (A < B) + (A < C))

        alpha = Symbol("alpha")
        for A, B in product(A_blades, B_blades):
            self.assertEquals((alpha * A) < B, alpha * (A < B))
            self.assertEquals((alpha * A) < B, A < (alpha * B))

        a = R.mv('a', 1, 'grade')
        for A_minus1, B in product(A_blades[:-1], B_blades):
            A = A_minus1 ^ a
            self.assertEquals(A < B, (A_minus1 ^ a) < B)
            self.assertEquals(A < B, A_minus1 < (a < B))


    def test_3_4(self):
        """
        The other contraction.
        """

        Ga.dual_mode("Iinv+")

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]
        B_blades = [R.mv('B', i, 'grade') for i in range(R.n + 1)]

        for A, B in product(A_blades, B_blades):
            self.assertEquals(B > A, ((-1) ** (A.pure_grade() * (B.pure_grade() - 1))) * (A < B))

        #for A, B in product(A_blades, B_blades):
        #    self.assertEquals((B > A).pure_grade(), B.pure_grade() - A.pure_grade())


    def test_3_5_2(self):
        """
        The inverse of a blade.
        """

        Ga.dual_mode("Iinv+")

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]

        for A in A_blades:
            self.assertEquals(A.inv(), ((-1) ** (A.pure_grade() * (A.pure_grade() - 1) / 2)) * (A / A.norm2()))

        for A in A_blades:
            self.assertEquals(A < A.inv(), 1)

        A = A_blades[1]
        self.assertEquals(A.inv(), A / A.norm2())


    def test_3_5_3(self):
        """
        Orthogonal complement and duality.
        """

        Ga.dual_mode("Iinv+")

        # some blades by grades for each space
        spaces = [([R.mv('A', i, 'grade') for i in range(R.n + 1)], R) for R in [Ga('e*1|2'), Ga('e*1|2|3'), Ga('e*1|2|3|4')]]
        
        # dualization
        for blades, R in spaces:
            for A in blades:
                self.assertEquals(A.dual(), A < R.I_inv())

        # dualization sign
        for blades, R in spaces:
            for A in blades:
                self.assertEquals(A.dual().dual(), ((-1) ** (R.n * (R.n - 1) / 2)) * A)

        # undualization
        for blades, R in spaces:
            for A in blades:
                self.assertEquals(A, A.dual() < R.I())


    def test_3_5_4(self):
        """
        The duality relationships.
        """

        Ga.dual_mode("Iinv+")

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]
        B_blades = [R.mv('B', i, 'grade') for i in range(R.n + 1)]

        for A, B in product(A_blades, B_blades):
            self.assertEquals((A ^ B).dual(), A < B.dual())

        for A, B in product(A_blades, B_blades):
            self.assertEquals((A < B).dual(), A ^ B.dual())


    def test_3_6(self):
        """
        Orthogonal projection of subspaces.
        """

        Ga.dual_mode("Iinv+")

        R = Ga('e*1|2|3')
        X_blades = [R.mv('X', i, 'grade') for i in range(R.n + 1)]
        B_blades = [R.mv('B', i, 'grade') for i in range(R.n + 1)]

        # projection of X on B
        def P(X, B):
            return (X < B.inv()) < B

        # a projection should be idempotent
        for X, B in product(X_blades, B_blades):
            self.assertEquals(P(X, B), P(P(X, B), B))

        # with the contraction
        for X, B in product(X_blades, B_blades):
            self.assertEquals(X < B, P(X, B) < B)
            
    
    def test_3_10_1_1(self):
        """
        Drills.
        """
        
        Ga.dual_mode("Iinv+")

        _R, e_1, e_2, e_3 = Ga.build('e*1|2|3', g='1 0 0, 0 1 0, 0 0 1')
        
        a = e_1 + e_2
        b = e_2 + e_3
        
        # a
        self.assertEquals(e_1 < a, (e_1 < e_1) + (e_1 < e_2))
        self.assertEquals(e_1 < e_1, e_1 | e_1)
        self.assertEquals(e_1 < e_1, 1)
        self.assertEquals(e_1 < e_2, e_1 | e_2)
        self.assertEquals(e_1 < e_2, 0)
        self.assertEquals(e_1 < a, 1)
        
        # b
        self.assertEquals(e_1 < (a ^ b), ((e_1 < a) ^ b) - (a ^ (e_1 < b)))        
        self.assertEquals((e_1 < a) ^ b, b)         # reusing drill a
        self.assertEquals(e_1 < b, (e_1 < e_2) + (e_1 < e_3))
        self.assertEquals(e_1 < b, 0)
        self.assertEquals(e_1 < (a ^ b), b)
        
        # c
        self.assertEquals((a ^ b) < e_1, a < (b < e_1))
        self.assertEquals((a ^ b) < e_1, a < ((e_2 < e_1) + (e_3 < e_1)))
        self.assertEquals((a ^ b) < e_1, 0)
        
        # d
        self.assertEquals(((2 * a) + b) < (a + b), ((2 * a) < a) + ((2 * a) < b) + (b < a) + (b < b))
        self.assertEquals(((2 * a) + b) < (a + b), ((2 * (a < a)) + (2 * (a < b)) + (b < a) + (b < b)))
        self.assertEquals(((2 * a) + b) < (a + b), 4 + 2 + 1 + 2)
        
        # e
        self.assertEquals(a < (e_1 ^ e_2 ^ e_3), a < ((e_1 ^ e_2) ^ e_3))
        self.assertEquals(a < (e_1 ^ e_2 ^ e_3), ((a < (e_1 ^ e_2)) ^ e_3) + ((e_1 ^ e_2) ^ (a < e_3)))
        self.assertEquals(a < (e_1 ^ e_2 ^ e_3), (((e_1 < (e_1 ^ e_2)) + (e_2 < (e_1 ^ e_2))) ^ e_3) + ((e_1 ^ e_2) ^ ((e_1 < e_3) + (e_2 < e_3))))        
        self.assertEquals(((e_1 < (e_1 ^ e_2)) + (e_2 < (e_1 ^ e_2))) ^ e_3, (((e_1 < e_1) ^ e_2) - (e_1 ^ (e_1 < e_2)) + ((e_2 < e_1) ^ e_2) - (e_1 ^ (e_2 < e_2))) ^ e_3)        
        self.assertEquals(((e_1 < (e_1 ^ e_2)) + (e_2 < (e_1 ^ e_2))) ^ e_3, (e_2 - e_1) ^ e_3)
        self.assertEquals(((e_1 < (e_1 ^ e_2)) + (e_2 < (e_1 ^ e_2))) ^ e_3, (e_2 ^ e_3) - (e_1 ^ e_3))        
        self.assertEquals(e_1 < e_3, 0)
        self.assertEquals(e_2 < e_3, 0)
        self.assertEquals(((e_1 ^ e_2) ^ ((e_1 < e_3) + (e_2 < e_3))), 0)        
        self.assertEquals(a < (e_1 ^ e_2 ^ e_3), (e_2 ^ e_3) - (e_1 ^ e_3))


if __name__ == '__main__':

    unittest.main()
