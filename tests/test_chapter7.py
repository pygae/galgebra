import unittest

from sympy import simplify, Symbol, pi, cos, sin, solve, sqrt, Rational, Mod
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

    def test7_9_1(self):
        """
        Drills.
        """
        Ga.dual_mode("Iinv+")

        GA, e_1, e_2, e_3 = Ga.build("e*1|2|3", g='1 0 0, 0 1 0, 0 0 1')

        # .1 Compute R1
        R1 = ((e_1 ^ e_2) * (-pi / 4)).exp()
        self.assertEquals(R1 * e_1 * R1.rev(), e_2)

        # .2 Compute R2
        R2 = ((e_3 ^ e_1) * (pi / 4)).exp()

        # .3 Compute R2 R1
        R2R1 = R2 * R1
        self.assertEquals(R2R1 * (e_1 ^ e_2) * R2R1.rev(), -e_2 ^ e_3)

        # .4 Compute the axis and angle of R2 R1
        a_1 = Symbol('a_1')
        a_2 = Symbol('a_2')
        a_3 = Symbol('a_3')
        a = GA.mv((a_1, a_2, a_3), 'vector')

        theta = Symbol('theta')
        A = a.dual()
        R = cos(theta / 2) + A * sin(theta / 2)

        system = (R2R1 - R).blade_coefs()
        unknowns = [a_1, a_2, a_3, theta]

        results = solve(system, unknowns, dict=True)
        self.assertEquals(a.subs(results[0]), (-e_1 - e_2 + e_3) / sqrt(3))
        self.assertEquals(Mod(theta.subs(results[0]), 2 * pi), Mod(Rational(2, 3) * pi, 2 * pi))
        self.assertEquals(a.subs(results[1]), -(-e_1 - e_2 + e_3) / sqrt(3))
        self.assertEquals(Mod(theta.subs(results[1]), 2 * pi), Mod(-Rational(2, 3) * pi, 2 * pi))

        # Check .4 against .3
        a = a.subs(results[0])
        theta = theta.subs(results[0])
        A = a.dual()
        R = cos(theta / 2) + A * sin(theta / 2)
        self.assertEquals(R * (e_1 ^ e_2) * R.rev(), -e_2 ^ e_3)

        # .5 Compute the product of rotors
        GA, e_1, e_2, e_3, e_4 = Ga.build("e*1|2|3|4", g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        R1 = ((e_1 ^ e_4) * (-pi / 4)).exp()
        R2 = ((e_2 ^ e_3) * (-pi / 4)).exp()
        R = R2 * R1

        self.assertEquals(R * (e_1 ^ e_2) * R.rev(), -e_3 ^ e_4)

        # .6 Reflect in a plane
        A = (e_1 + e_2) ^ e_3
        B = e_1 ^ e_4
        b = B.dual()

        def hat(k, M):
            return ((-1) ** k) * M

        self.assertEquals(b * hat(2, A) * b.inv(), (-e_1 + e_2) ^ e_3)

        # .7 Reflect the dual plane reflector e1 in the plane e1 ^ e3
        GA, e_1, e_2, e_3 = Ga.build("e*1|2|3", g='1 0 0, 0 1 0, 0 0 1')

        a = e_1
        B = e_1 ^ e_3
        b = B.dual()

        self.assertEquals(-b * a * b.inv(), e_1)

        # Reset
        Ga.dual_mode()


if __name__ == '__main__':

    unittest.main()
