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


    def test_3_1_2(self):
        """
        Definition of the scalar product.
        """

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]
        B_blades = [R.mv('B', i, 'grade') for i in range(R.n + 1)]

        # GAlgebra doesn't define any scalar product but rely on the geometric product instead
        for A, B in product(A_blades, B_blades):
            A_grades = R.grade_decomposition(A).keys()
            B_grades = R.grade_decomposition(B).keys()
            self.assertEquals(len(A_grades), 1)
            self.assertEquals(len(B_grades), 1)
            if A_grades[0] == B_grades[0]:
                self.assertTrue((A * B).scalar() != 0)
            else:
                self.assertTrue((A * B).scalar() == 0)


    def test_3_2_2(self):
        """
        Computing the contraction explicitly.
        """

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
            B_grades = R.grade_decomposition(B).keys()
            self.assertEquals(len(B_grades), 1)
            self.assertEquals(B < A, 0 if B_grades[0] else A * B)

        # vectors
        A = A_blades[1]
        B = B_blades[1]
        self.assertEquals(A < B, A | B)

        # vector and the outer product of 2 blades of various grades (scalars, vectors, 2-vectors...)
        A = A_blades[1]
        for B, C in product(B_blades, C_blades):
            B_grades = R.grade_decomposition(B).keys()
            self.assertEquals(len(B_grades), 1)
            self.assertEquals(A < (B ^ C), ((A < B) ^ C) + (-1)**B_grades[0] * (B ^ (A < C)))

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

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]
        B_blades = [R.mv('B', i, 'grade') for i in range(R.n + 1)]

        for A, B in product(A_blades, B_blades):
            A_grades = R.grade_decomposition(A).keys()
            B_grades = R.grade_decomposition(B).keys()
            self.assertEquals(len(A_grades), 1)
            self.assertEquals(len(B_grades), 1)
            self.assertEquals(B > A, ((-1) ** (A_grades[0] * (B_grades[0] - 1))) * (A < B))

        for A, B in product(A_blades, B_blades):
            C = B > A
            A_grades = R.grade_decomposition(A).keys()
            B_grades = R.grade_decomposition(B).keys()
            C_grades = R.grade_decomposition(C).keys()
            self.assertEquals(len(A_grades), 1)
            self.assertEquals(len(B_grades), 1)
            self.assertEquals(len(C_grades), 1)
            self.assertTrue(C == 0 or C_grades[0] == B_grades[0] - A_grades[0])


    def test_3_5_2(self):
        """
        The inverse of a blade.
        """

        R = Ga('e*1|2|3')
        A_blades = [R.mv('A', i, 'grade') for i in range(R.n + 1)]

        for A in A_blades:
            A_grades = R.grade_decomposition(A).keys()
            self.assertEquals(len(A_grades), 1)
            self.assertEquals(A.inv(), ((-1) ** (A_grades[0] * (A_grades[0] - 1) / 2)) * (A / A.norm2()))

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


    def test3_7_2(self):
        """
        The cross product incorporated.
        """

        Ga.dual_mode("Iinv+")

        R, e_1, e_2, e_3 = Ga.build('e*1|2|3')
        a = R.mv('a', 'vector')
        b = R.mv('b', 'vector')

        A = a < R.I()
        B = b < R.I()

        # Cross product definition
        self.assertEquals(R.I_inv(), -R.I())
        self.assertEquals(cross(a, b), (a ^ b).dual())

        # About velocities
        self.assertEquals(cross(a, b), (b ^ a) < R.I())
        self.assertEquals(cross(a, b), (a ^ b) < R.I_inv())
        self.assertEquals(cross(a, b), -(b ^ a).dual())
        self.assertEquals(cross(a, b), -b < a.dual())
        self.assertEquals(cross(a, b), b < A)

        # Intersecting planes
        self.assertEquals(cross(a, b), ((A < R.I_inv()) ^ (B < R.I_inv())) < R.I_inv())
        self.assertEquals(cross(a, b), (B < R.I_inv()) < ((A < R.I_inv()) < R.I()))
        self.assertEquals(cross(a, b), (B < R.I_inv()) < A)
        self.assertEquals(cross(a, b), B.dual() < A)

    
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
        
        R, e_1, e_2, e_3, e_4 = Ga.build('e*1|2|3|4', g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')
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


    def test3_10_2_1(self):
        """
        In 2-D Euclidean space R(2,0) with orthogonal basis {e1, e2}, let us determine the value of the contraction
        e1 < (e1 ^ e2) by means of its implicit definition with A = e1 and B = e1 ^ e2. Let X range over the basis of
        the blades {1, e1, e2, e1 ^ e2}. This produces four equations, each of which gives you information on the
        coefficient of the corresponding basis element in the final result.
        Show that e1 < (e1 ^ e2) = (0)1 + 0(e1) + 1(e2) + 0(e1 ^ e2).
        """

        R, e_1, e_2 = Ga.build('e*1|2', g='1 0, 0 1')

        A = e_1
        B = e_1 ^ e_2

        X = R.mv(1)
        self.assertEquals(((X ^ A) * B).scalar(), (X * (A < B)).scalar())
        self.assertEquals(((X ^ A) * B).scalar(), 0)

        X = e_1
        self.assertEquals(((X ^ A) * B).scalar(), (X * (A < B)).scalar())
        self.assertEquals(((X ^ A) * B).scalar(), 0)

        X = e_2
        self.assertEquals(((X ^ A) * B).scalar(), (X * (A < B)).scalar())
        self.assertEquals(((X ^ A) * B).scalar(), 1)

        X = e_1 ^ e_2
        self.assertEquals(((X ^ A) * B).scalar(), (X * (A < B)).scalar())
        self.assertEquals(((X ^ A) * B).scalar(), 0)


    def test3_10_2_2(self):
        """
        Change the metric such that e2 . e2 == 0. Show that you can't determine the coefficient of e2 in the value
        of e1 < (e1 ^ e2) like the previous exercise. The use the explicit definition of the contraction to show the
        contraction is still well defined, and equal to e1 < (e1 ^ e2) == e2.
        """

        R, e_1, e_2 = Ga.build('e*1|2', g='1 0, 0 0')

        # We can't use the scalar product anymore (because of the metric)
        A = e_1
        B = e_1 ^ e_2
        X = e_2
        self.assertEquals(((X ^ A) * B).scalar(), (X * (A < B)).scalar())
        self.assertEquals(((X ^ A) * B).scalar(), 0)    # Which is false

        e_1_grades = R.grade_decomposition(e_1).keys()
        self.assertEquals(len(e_1_grades), 1)
        self.assertEquals(e_1_grades[0], 1)

        # We use the explicit definition of the contraction instead
        self.assertEquals(e_1 < (e_1 ^ e_2), ((e_1 < e_1) ^ e_2) + ((-1) ** e_1_grades[0]) * (e_1 ^ (e_1 < e_2)))
        # We can't use the definition a < b = a . b using GAlgebra so we solved it ourselves...
        self.assertEquals(e_1 < (e_1 ^ e_2), (R.mv(1) ^ e_2) + ((-1) ** e_1_grades[0]) * (e_1 ^ R.mv(0)))
        self.assertEquals(e_1 < (e_1 ^ e_2), e_2)       # Which is true


    def test3_10_2_3(self):
        """
        Derive the following dualities for right contraction :
        C > (B ^ A) = (C > B) > A and C > (B > A) = (C > B) ^ A when A included in C.
        """

        R, e_1, e_2, e_3 = Ga.build('e*1|2|3')

        A_blades = [R.mv('X', k, 'grade') for k in range(R.n + 1)]
        B_blades = [R.mv('B', l, 'grade') for l in range(R.n + 1)]
        C_blades = [R.mv('B', m, 'grade') for m in range(R.n + 1)]

        for A, B, C in product(A_blades, B_blades, C_blades):
            M = C > (B ^ A)
            self.assertEquals(M, ((A.rev() ^ B.rev()) < C.rev()).rev())
            self.assertEquals(M, (A.rev() < (B.rev() < C.rev())).rev())
            self.assertEquals(M, (B.rev() < C.rev()).rev() > A)
            self.assertEquals(M, (C > B) > A)

            M = C > (B > A)
            # TODO

    def test3_10_2_11(self):
        """
        Derive the notorious bac-cab formula for the cross product (a x (b x c) = b (a . c) - c (a . b)),
        directly from its definition (3.28). What is the corresponding formula using ^ and < ?
        """
        Ga.dual_mode("Iinv+")

        R, e_1, e_2, e_3 = Ga.build('e*1|2|3')
        a = R.mv('a', 'vector')
        b = R.mv('b', 'vector')
        c = R.mv('c', 'vector')

        xx = cross(a, cross(b, c))
        self.assertEquals(xx, (a ^ (b ^ c).dual()).dual())
        self.assertEquals(xx, a < (b ^ c).dual().dual())
        self.assertEquals(xx, -a < (b ^ c))
        self.assertEquals(xx, (-a < b) ^ c - b ^ (-a < c))
        self.assertEquals(xx, b ^ (a < c) - c ^ (a < b))
        self.assertTrue((a < c).is_scalar())
        self.assertTrue((a < b).is_scalar())
        self.assertEquals(xx, b * (a < c) - c * (a < b))


if __name__ == '__main__':

    unittest.main()
