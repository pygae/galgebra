from .test_utils import TestCase

from sympy import Symbol, pi, cos, sin, solve, sqrt, Rational, Mod
from galgebra.ga import Ga


class TestChapter7(TestCase):

    def test7_9_1(self):
        """
        Drills.
        """
        Ga.dual_mode("Iinv+")

        GA, e_1, e_2, e_3 = Ga.build("e*1|2|3", g='1 0 0, 0 1 0, 0 0 1')

        # .1 Compute R1
        R1 = ((e_1 ^ e_2) * (-pi / 4)).exp()
        self.assertEqual(R1 * e_1 * R1.rev(), e_2)

        # .2 Compute R2
        R2 = ((e_3 ^ e_1) * (pi / 4)).exp()

        # .3 Compute R2 R1
        R2R1 = R2 * R1
        self.assertEqual(R2R1 * (e_1 ^ e_2) * R2R1.rev(), -e_2 ^ e_3)

        # .4 Compute the axis and angle of R2 R1
        a_1 = Symbol('a_1', real=True)
        a_2 = Symbol('a_2', real=True)
        a_3 = Symbol('a_3', real=True)
        a = GA.mv((a_1, a_2, a_3), 'vector')

        theta = Symbol('theta', real=True)
        A = a.dual()
        R = cos(theta / 2) + A * sin(theta / 2)

        system = (R2R1 - R).blade_coefs()
        unknowns = [a_1, a_2, a_3, theta]

        results = solve(system, unknowns, dict=True)
        self.assertEqual(a.subs(results[0]), (-e_1 - e_2 + e_3) / sqrt(3))
        self.assertEqual(Mod(theta.subs(results[0]), 2 * pi), Mod(Rational(2, 3) * pi, 2 * pi))
        self.assertEqual(a.subs(results[1]), -(-e_1 - e_2 + e_3) / sqrt(3))
        self.assertEqual(Mod(theta.subs(results[1]), 2 * pi), Mod(-Rational(2, 3) * pi, 2 * pi))

        # Check .4 against .3
        a = a.subs(results[0])
        theta = theta.subs(results[0])
        A = a.dual()
        R = cos(theta / 2) + A * sin(theta / 2)
        self.assertEqual(R * (e_1 ^ e_2) * R.rev(), -e_2 ^ e_3)

        # .5 Compute the product of rotors
        GA, e_1, e_2, e_3, e_4 = Ga.build("e*1|2|3|4", g='1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        R1 = ((e_1 ^ e_4) * (-pi / 4)).exp()
        R2 = ((e_2 ^ e_3) * (-pi / 4)).exp()
        R = R2 * R1

        self.assertEqual(R * (e_1 ^ e_2) * R.rev(), -e_3 ^ e_4)

        # .6 Reflect in a plane
        A = (e_1 + e_2) ^ e_3
        B = e_1 ^ e_4
        b = B.dual()

        def hat(k, M):
            return ((-1) ** k) * M

        self.assertEqual(b * hat(2, A) * b.inv(), (-e_1 + e_2) ^ e_3)

        # .7 Reflect the dual plane reflector e1 in the plane e1 ^ e3
        GA, e_1, e_2, e_3 = Ga.build("e*1|2|3", g='1 0 0, 0 1 0, 0 0 1')

        a = e_1
        B = e_1 ^ e_3
        b = B.dual()

        self.assertEqual(-b * a * b.inv(), e_1)

        # Reset
        Ga.dual_mode()
