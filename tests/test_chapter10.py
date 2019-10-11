import unittest

from sympy import simplify, sqrt, Rational, Symbol
from ga import Ga
from mv import Mv


class TestChapter7(unittest.TestCase):

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

    def test_11_12_1(self):
        """
        Compute the 2-blades corresponding to the lines gives by the data below. Which of
        the lines are the same, considered as weighted oriented elements of geometry, which
        are the same as offset subspaces ?
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        p = e_0 + e_1
        q = e_0 + e_2
        d = e_2 - e_1
        self.assertEquals(p ^ q, p ^ d)
        self.assertEquals(p ^ q, q ^ d)

        r = e_0 + 2 * (e_2 - e_1)
        e = 2 * (e_2 - e_1)
        s = e_0 + 3 * (e_2 - e_1)
        t = 2 * (e_0 + e_2)
        self.assertEquals(2 * (p ^ q), p ^ e)
        self.assertEquals(2 * (p ^ q), p ^ t)

    def test11_12_2_1(self):
        """
        Let an orthonormal coordinate system be given in 3-dimensional Euclidean space.
        Compute the support vector of the line with direction u = e1 + 2 e2 - e3, through the point p = e1 - 3 e2.
        What is the distance of the line to the origin ?
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        L = (e_1 + 2 * e_2 - e_3) ^ (e_0 + e_1 - 3 * e_2)
        d = (e_0.inv() < (e_0 ^ L)) / (e_0.inv() < L)
        self.assertEquals(d.norm(), sqrt(Rational(35, 6)))

    def test11_12_2_2(self):
        """
        Convert the line of the previous exercise into a parametric equation x = p + t * u; express t as a function of x
        for a point x on the line.
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        u = e_1 + 2 * e_2 - e_3
        p = e_0 + e_1 - 3 * e_2
        L = u ^ p

        t = Symbol('t')
        x_1 = Symbol('x_1')
        x_2 = Symbol('x_2')
        x_3 = Symbol('x_3')
        x = e_0 + x_1 * e_1 + x_2 * e_2 + x_3 * e_3

        # x(t)
        x_t = (e_0 + e_1 - 3 * e_2) + t * (e_1 + 2 * e_2 - e_3)

        # t(x)
        t_x = (x - (e_0 + e_1 - 3 * e_2)) / (e_1 + 2 * e_2 - e_3)

        for i in range(11):
            t_value = Rational(i, 10)
            x_value = x_t.subs({t: t_value})
            self.assertEquals(x_value ^ L, 0)
            self.assertEquals(t_x.subs(zip([x_1, x_2, x_3], x_value.blade_coefs([e_1, e_2, e_3]))), t_value)
