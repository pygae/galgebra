import unittest

from sympy import S, simplify
from galgebra.ga import Ga
from galgebra.mv import Mv


def com(A, B):
    """
    I like free functions...
    """
    return Ga.com(A, B)


class TestCase(unittest.TestCase):

    def assertEquals(self, first, second):
        """
        Compare two expressions are equals.
        """
        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        diff = simplify(first - second)

        self.assertTrue(diff == 0, "\n%s\n==\n%s\n%s" % (first, second, diff))

    def assertProjEquals(self, X, Y):
        """
        Compare two points, two planes or two lines up to a scalar.
        """
        assert isinstance(X, Mv)
        assert isinstance(Y, Mv)

        X /= self.norm(X)
        Y /= self.norm(Y)

        # We can't easily retrieve the sign, so we test both
        diff = simplify(X.obj - Y.obj)
        if diff != S.Zero:
            diff = simplify(X.obj + Y.obj)

        self.assertTrue(diff == S.Zero, "\n%s\n==\n%s" % (X, Y))

    def assertNotEquals(self, first, second):
        """
        Compare two expressions are not equals.
        """
        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        diff = simplify(first - second)

        self.assertTrue(diff != 0, "\n%s\n!=\n%s\n%s" % (first, second, diff))
