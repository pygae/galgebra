from .test_utils import TestCase
import importlib

from galgebra.ga import Ga
from galgebra.mv import J, Jinv
from galgebra.generator import format_geometric_algebra, flatten, expand


class TestGenerator(TestCase):

    def setUp(self):

        Ga.dual_mode("Iinv+")
        GA = Ga('e*1|2|3|4', g=[0, 1, 1, 1])
        GA.build_cobases()
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

    def test_and(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(Jinv(J(X) ^ J(Y)), expand(GA, flatten(flat_GA, X) & flatten(flat_GA, Y)))

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

    def test_meet_and_join(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')
        Y = GA.mv('y', 'mv')

        self.assertEquals(Jinv(J(X) ^ J(Y)), expand(GA, flatten(flat_GA, X).meet(flatten(flat_GA, Y))))
        self.assertEquals(X ^ Y, expand(GA, flatten(flat_GA, X).join(flatten(flat_GA, Y))))

    def test_rev(self):

        GA = self.GA
        flat_GA = self.flat_GA
        X = GA.mv('x', 'mv')

        self.assertEquals(X.rev(), expand(GA, flatten(flat_GA, X).rev()))
