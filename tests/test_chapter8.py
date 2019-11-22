import unittest

from sympy import simplify, Symbol, S

from galgebra.ga import Ga
from galgebra.mv import Mv, com


class TestChapter8(unittest.TestCase):

    def assertEquals(self, first, second, msg=None):
        """
        Compare two expressions are equals.
        """

        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        diff = simplify(first - second)

        self.assertTrue(diff == 0, "\n%s\n==\n%s\n%s" % (first, second, diff))

    def test8_2_1(self):
        """
        The commutator product.
        """
        GA = Ga("e*1|2|3", g=[1, 1, 1])

        A = GA.mv('A', 'mv')
        B = GA.mv('B', 'mv')
        C = GA.mv('C', 'mv')

        self.assertEquals(com(com(A, B), C) - com(A, com(B, C)), com(B, com(C, A)))
        self.assertEquals(com(com(A, B), C) + com(com(C, A), B) + com(com(B, C), A), S.Zero)

        B = GA.mv('B', 2, 'blade')
        X = GA.mv('A', 'mv')
        self.assertEquals(B.pure_grade(), 2)
        self.assertEquals(com(X, B).pure_grade(), -2)   # Not pure

        E = com(X, B).rev()
        self.assertEquals(E, (B.rev() * X.rev() - X.rev() * B.rev()) / 2)
        self.assertEquals(E, (X.rev() * B - B * X.rev()) / 2)
        self.assertEquals(E, com(X.rev(), B))

        alpha = GA.mv('alpha', 'scalar')
        self.assertEquals(alpha * X, alpha ^ X)

        a = GA.mv('a', 'vector')
        self.assertEquals(a * X, (a < X) + (a ^ X))

        A = GA.mv('A', 2, 'grade')
        self.assertEquals(A * X, (A < X) + com(A, X) + (A ^ X))
