import unittest
import importlib

from sympy import simplify

from ga import Ga
from mv import Mv
from generator import format_geometric_algebra, flatten, expand


class TestGenerator(unittest.TestCase):

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

    def setUp(self):

        Ga.dual_mode("Iinv+")
        GA = Ga('e*1|2|3', g=[1, 1, 1])
        with open('flat_ga.py', 'wt') as flat_GA_file:
            flat_GA_file.write(format_geometric_algebra(GA))
        self.GA = GA

        flat_GA = importlib.import_module('flat_ga')
        self.flat_GA = flat_GA

    def test_add(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(X + Y, expand(GA, flatten(flat_GA, X) + flatten(flat_GA, Y)))

    def test_sub(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(X - Y, expand(GA, flatten(flat_GA, X) - flatten(flat_GA, Y)))

    def test_mul(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(X * Y, expand(GA, flatten(flat_GA, X) * flatten(flat_GA, Y)))

    def test_xor(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(X ^ Y, expand(GA, flatten(flat_GA, X) ^ flatten(flat_GA, Y)))

    def test_lshift(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(X << Y, expand(GA, flatten(flat_GA, X) << flatten(flat_GA, Y)))

    def test_rshift(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(X >> Y, expand(GA, flatten(flat_GA, X) >> flatten(flat_GA, Y)))

    def test_dual(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')

        self.assertEquals(X.dual(), expand(GA, flatten(flat_GA, X).dual()))

    def test_rev(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')

        self.assertEquals(X.rev(), expand(GA, flatten(flat_GA, X).rev()))
