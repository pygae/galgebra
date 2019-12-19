from .test_utils import TestCase, com

from sympy import S
from galgebra.ga import Ga


class TestChapter8(TestCase):

    def test8_2_1(self):
        """
        The commutator product.
        """
        GA = Ga("e*1|2|3", g=[1, 1, 1])

        A = GA.mv('A', 'mv')
        B = GA.mv('B', 'mv')
        C = GA.mv('C', 'mv')

        self.assertEqual(com(com(A, B), C) - com(A, com(B, C)), com(B, com(C, A)))
        self.assertEqual(com(com(A, B), C) + com(com(C, A), B) + com(com(B, C), A), S.Zero)

        B = GA.mv('B', 'blade', 2)
        X = GA.mv('A', 'mv')
        self.assertEqual(B.pure_grade(), 2)
        self.assertEqual(com(X, B).pure_grade(), -2)  # Not pure

        E = com(X, B).rev()
        self.assertEqual(E, (B.rev() * X.rev() - X.rev() * B.rev()) / 2)
        self.assertEqual(E, (X.rev() * B - B * X.rev()) / 2)
        self.assertEqual(E, com(X.rev(), B))

        alpha = GA.mv('alpha', 'scalar')
        self.assertEqual(alpha * X, alpha ^ X)

        a = GA.mv('a', 'vector')
        self.assertEqual(a * X, (a < X) + (a ^ X))

        A = GA.mv('A', 'grade', 2)
        self.assertEqual(A * X, (A < X) + com(A, X) + (A ^ X))
