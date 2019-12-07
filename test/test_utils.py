import unittest

from sympy import expand, S, simplify
from galgebra.ga import Ga
from galgebra.mv import Mv
from galgebra.metric import Simp


def com(A, B):
    """
    I like free functions...
    """
    return Ga.com(A, B)


class TestCase(unittest.TestCase):

    def assertEqual(self, first, second):
        """
        Compare two expressions are equals.
        """
        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        # We need to help sympy a little...
        first = Simp.apply(expand(first))
        second = Simp.apply(expand(second))

        # Check
        diff = simplify(first - second)

        self.assertTrue(diff == 0, "\n%s\n==\n%s\n%s" % (first, second, diff))

    def assertProjEqual(self, first, second):
        """
        Compare two points, two planes or two lines up to a scalar.
        """
        assert isinstance(first, Mv)
        assert isinstance(second, Mv)

        # TODO: this should use Mv methods and not the derived test case methods...
        first /= self.norm(first)
        second /= self.norm(second)

        # We need to help sympy a little...
        X = Simp.apply(expand(first.obj))
        Y = Simp.apply(expand(second.obj))

        # We can't easily retrieve the sign, so we test both
        diff = simplify(X.obj - Y.obj)
        if diff != S.Zero:
            diff = simplify(X.obj + Y.obj)

        self.assertTrue(diff == S.Zero, "\n%s\n==\n%s" % (X, Y))

    def assertNotEqual(self, first, second):
        """
        Compare two expressions are not equals.
        """
        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        # We need to help sympy a little...
        first = Simp.apply(expand(first))
        second = Simp.apply(expand(second))

        # Check
        diff = simplify(first - second)

        self.assertTrue(diff != 0, "\n%s\n!=\n%s\n%s" % (first, second, diff))
