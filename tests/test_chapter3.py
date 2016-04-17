import unittest

from itertools import product
from sympy import simplify, Symbol, cos, sin, sqrt
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
        Let a = e_1 + e_2 and b = e_2 + e_3 in  a 3D Euclidean space R(3,0) with an orthonormal basis {e_1, e_2, e_3}.
        Compute the following expressions, giving the results relative to the basis
        {1, e_1, e_2, e_3, e_1 ^ e_2, e_2 ^ e_3, e_3 ^ e_1, e_1 ^ e_2 ^ e_3}.
        """
        
        Ga.dual_mode("Iinv+")

        R, e_1, e_2, e_3 = Ga.build('e*1|2|3', g='1 0 0, 0 1 0, 0 0 1')
        
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
        self.assertEquals(e_1 < (a ^ b), e_2 + e_3)
        
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
        self.assertEquals(a < (e_1 ^ e_2 ^ e_3), (e_2 ^ e_3) + (e_3 ^ e_1))
        
        # f
        self.assertEquals(a < R.I_inv(), a < (-e_1 ^ e_2 ^ e_3))        # equals because we fixed the metric
        self.assertEquals(a < R.I_inv(), a < ((e_1 ^ e_2) ^ (-e_3)))    # almost last drill
        self.assertEquals(a < R.I_inv(), -(e_2 ^ e_3) - (e_3 ^ e_1))
        
        # g
        self.assertEquals((a ^ b) < R.I_inv(), -((b ^ a) < R.I_inv()))
        self.assertEquals((a ^ b) < R.I_inv(), -(b < (a < R.I_inv())))
        self.assertEquals((a ^ b) < R.I_inv(), -((e_2 + e_3) < (-(e_2 ^ e_3) + (e_1 ^ e_3))))        
        self.assertEquals((a ^ b) < R.I_inv(), ((e_2 + e_3) < (e_2 ^ e_3)) - ((e_2 + e_3) < (e_1 ^ e_3)))
        self.assertEquals((a ^ b) < R.I_inv(), ((e_2 < (e_2 ^ e_3)) + (e_3 < (e_2 ^ e_3)) - (e_2 < (e_1 ^ e_3)) - (e_3 < (e_1 ^ e_3))))
        self.assertEquals(e_2 < (e_2 ^ e_3), e_3)
        self.assertEquals(e_3 < (e_2 ^ e_3), -e_2)
        self.assertEquals(e_2 < (e_1 ^ e_3), 0)
        self.assertEquals(e_3 < (e_1 ^ e_3), -e_1)
        self.assertEquals((a ^ b) < R.I_inv(), e_3 - e_2 + e_1)
        
        # h
        self.assertEquals(a < (b < R.I_inv()), (a ^ b) < R.I_inv())
        self.assertEquals(a < (b < R.I_inv()), e_3 - e_2 + e_1)
    
    
    def test_3_10_1_2(self):
        """
        Compute the cosine of the angle between the following subspaces given on an orthonormal basis of a Euclidean space.        
        """
                
        Ga.dual_mode("Iinv+")
        
        _R, e_1, e_2, e_3, e_4 = Ga.build('e*1|2|3|4', g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')                
        alpha = Symbol('alpha', real=True)
        theta = Symbol('theta', real=True)
        
        # e_1 and alpha * e_1                
        A = e_1
        B = alpha * e_1
        cosine = (A | B.rev()) / (A.norm() * B.norm())        
        self.assertEquals(A | B.rev(), e_1 | (alpha * e_1))
        self.assertEquals(A.norm() * B.norm(), alpha)        
        self.assertEquals(cosine, (e_1 | (alpha * e_1)) / alpha)
        self.assertEquals(cosine, (alpha * (e_1 | e_1)) / alpha)
        self.assertEquals(cosine, 1)
        
        # (e_1 + e_2) ^ e_3 and e_1 ^ e_3        
        A = (e_1 + e_2) ^ e_3
        B = e_1 ^ e_3
        num = A | B.rev()
        den = A.norm() * B.norm()
        cosine = num / den          
        self.assertEquals(num, ((e_1 ^ e_3) + (e_2 ^ e_3)) | (e_3 ^ e_1))
        self.assertEquals(num, ((e_1 ^ e_3) | (e_3 ^ e_1)) + ((e_2 ^ e_3) | (e_3 ^ e_1)))
        self.assertEquals(((e_1 ^ e_3) | (e_3 ^ e_1)), 1)
        self.assertEquals(((e_2 ^ e_3) | (e_3 ^ e_1)), 0)
        self.assertEquals(num, 1)        
        self.assertEquals(den, ((e_1 + e_2) ^ e_3).norm() * (e_1 ^ e_3).norm())        
        self.assertEquals(den, ((e_1 ^ e_3) + (e_2 ^ e_3)).norm())        
        den2 = ((e_1 ^ e_3) + (e_2 ^ e_3)).norm2()        
        self.assertEquals(den2, ((e_1 ^ e_3) + (e_2 ^ e_3)) | ((e_1 ^ e_3) + (e_2 ^ e_3)).rev())        
        self.assertEquals(den2, ((e_1 ^ e_3) + (e_2 ^ e_3)) | ((e_3 ^ e_1) + (e_3 ^ e_2)))        
        self.assertEquals(den2, ((e_1 ^ e_3) | (e_3 ^ e_1)) + ((e_1 ^ e_3) | (e_3 ^ e_2)) + ((e_2 ^ e_3) | (e_3 ^ e_1)) + ((e_2 ^ e_3) | (e_3 ^ e_2)))
        self.assertEquals(den2, 1 + 0 + 0 + 1)        
        self.assertEquals(cosine, 1 / sqrt(2))
        
        # (cos(theta) * e_1 + sin(theta) * e_2) ^ e_3 and e_2 ^ e_3        
        A = (cos(theta) * e_1 + sin(theta) * e_2) ^ e_3
        B = e_2 ^ e_3
        num = A | B.rev()
        den = A.norm() * B.norm()
        cosine = num / den        
        self.assertEquals(num, (cos(theta) * (e_1 ^ e_3) + sin(theta) * (e_2 ^ e_3)) | (e_3 ^ e_2))
        self.assertEquals(num, ((cos(theta) * (e_1 ^ e_3)) | (e_3 ^ e_2)) + ((sin(theta) * (e_2 ^ e_3)) | (e_3 ^ e_2)))
        self.assertEquals((cos(theta) * (e_1 ^ e_3)) | (e_3 ^ e_2), 0)
        self.assertEquals((sin(theta) * (e_2 ^ e_3)) | (e_3 ^ e_2), sin(theta))
        self.assertEquals(num, sin(theta))
        self.assertEquals(den, ((cos(theta) * e_1 + sin(theta) * e_2) ^ e_3).norm() * (e_2 ^ e_3).norm())
        self.assertEquals(den, ((cos(theta) * e_1 + sin(theta) * e_2) ^ e_3).norm())        
        den2 = ((cos(theta) * e_1 + sin(theta) * e_2) ^ e_3).norm2()
        self.assertEquals(den2, ((cos(theta) * e_1 + sin(theta) * e_2) ^ e_3) | ((cos(theta) * e_1 + sin(theta) * e_2) ^ e_3).rev())
        self.assertEquals(den2, (cos(theta) * (e_1 ^ e_3) + sin(theta) * (e_2 ^ e_3)) | (cos(theta) * (e_1 ^ e_3) + sin(theta) * (e_2 ^ e_3)).rev())
        self.assertEquals(den2, (cos(theta) * (e_1 ^ e_3) + sin(theta) * (e_2 ^ e_3)) | (cos(theta) * (e_3 ^ e_1) + sin(theta) * (e_3 ^ e_2)))
        self.assertEquals(den2, ((cos(theta) * (e_1 ^ e_3)) | (cos(theta) * (e_3 ^ e_1))) + ((cos(theta) * (e_1 ^ e_3)) | (sin(theta) * (e_3 ^ e_2))) + ((sin(theta) * (e_2 ^ e_3)) | (cos(theta) * (e_3 ^ e_1))) + ((sin(theta) * (e_2 ^ e_3)) | (sin(theta) * (e_3 ^ e_2))))
        self.assertEquals(den2, cos(theta) * cos(theta) + sin(theta) * sin(theta))
        self.assertEquals(den2, 1)        
        self.assertEquals(cosine, sin(theta))
                
        # e_1 ^ e_2 and e_3 ^ e_4        
        A = e_1 ^ e_2
        B = e_3 ^ e_4
        num = A | B.rev()
        den = A.norm() * B.norm()
        cosine = num / den
        self.assertEquals(num, (e_1 ^ e_2) | (e_4 ^ e_3))
        self.assertEquals(num, 0)        
        self.assertEquals(cosine, 0)


if __name__ == '__main__':

    unittest.main()
