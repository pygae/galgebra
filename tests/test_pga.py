import unittest

from sympy import simplify, Rational, symbols, Symbol, S

from ga import Ga
from mv import Mv, J, Jinv


class TestChapter11(unittest.TestCase):

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
        """
        Setup 3D Projective Geometric Algebra aka PGA.
        """
        PGA, e_0, e_1, e_2, e_3 = Ga.build('e*0|1|2|3', g=[0, 1, 1, 1])

        # TODO: move this somewhere useful...
        PGA.build_cobases()

        self.PGA = PGA
        self.e_0 = e_0
        self.e_1 = e_1
        self.e_2 = e_2
        self.e_3 = e_3
        self.e_032 = e_0 ^ e_3 ^ e_2
        self.e_013 = e_0 ^ e_1 ^ e_3
        self.e_021 = e_0 ^ e_2 ^ e_1
        self.e_123 = e_1 ^ e_2 ^ e_3
        self.e_0123 = e_0 ^ e_1 ^ e_2 ^ e_3

    def homogenize(self, P):
        """
        For testing equality.
        """
        return P / (P.blade_coefs([self.e_123])[0])

    def point(self, x, y, z):
        """
        Make a point.
        """
        return x * self.e_032 + y * self.e_013 + z * self.e_021 + self.e_123

    def direction(self, x, y, z):
        """
        Make an ideal point (direction).
        """
        return x * self.e_032 + y * self.e_013 + z * self.e_021

    def plane(self, a, b, c, d):
        """
        Make a plane.
        """
        return a * self.e_1 + b * self.e_2 + c * self.e_3 + d * self.e_0

    def test_J(self):

        PGA = self.PGA

        for k in range(PGA.n + 1):
            X = PGA.mv('x', k, 'grade')
            self.assertEquals(X, Jinv(J(X)))

        X = PGA.mv('x', 'mv')
        self.assertEquals(X, Jinv(J(X)))

    def test_geometry_incidence_planes_meet_into_points(self):
        """
        Planes meet into points.
        """
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        p1 = self.plane(1, 0, 0, -x)
        p2 = self.plane(0, 1, 0, -y)
        p3 = self.plane(0, 0, 1, -z)

        P = self.homogenize(p1 ^ p2 ^ p3)
        self.assertEquals(P, self.point(x, y, z))
        self.assertEquals(P, self.homogenize(p1 ^ p3 ^ p2))
        self.assertEquals(P, self.homogenize(p2 ^ p1 ^ p3))
        self.assertEquals(P, self.homogenize(p2 ^ p3 ^ p1))
        self.assertEquals(P, self.homogenize(p3 ^ p1 ^ p2))
        self.assertEquals(P, self.homogenize(p3 ^ p2 ^ p1))

    def test_geometry_incidence_planes_meet_into_points_2(self):
        """
        Planes meet into points.
        """
        P1 = self.point(*symbols('x1 y1 z1'))
        P2 = self.point(*symbols('x2 y2 z2'))
        P3 = self.point(*symbols('x3 y3 z3'))
        P4 = self.point(*symbols('x4 y4 z4'))
        p123 = Jinv(J(P1) ^ J(P2) ^ J(P3))
        p124 = Jinv(J(P1) ^ J(P2) ^ J(P4))
        p234 = Jinv(J(P2) ^ J(P3) ^ J(P4))
        p314 = Jinv(J(P3) ^ J(P1) ^ J(P4))

        # TODO: find a way to meet planes faster...
        self.assertEquals(self.homogenize(p123 ^ p124 ^ p314), P1)
        self.assertEquals(self.homogenize(p123 ^ p124 ^ p234), P2)
        self.assertEquals(self.homogenize(p123 ^ p234 ^ p314), P3)
        self.assertEquals(self.homogenize(p124 ^ p234 ^ p314), P4)

    def test_geometry_incidence_points_join_into_planes(self):
        """
        Points join into planes.
        """
        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        x2, y2, z2 = symbols('x2 y2 z2')
        P2 = self.point(x2, y2, z2)

        x3, y3, z3 = symbols('x3 y3 z3')
        P3 = self.point(x3, y3, z3)

        pp = Jinv(J(P1) ^ J(P2) ^ J(P3))
        pm = Jinv(J(P1) ^ J(P3) ^ J(P2))
        self.assertEquals(pp, -pm)

        px = Symbol('px')
        py = Symbol('py')
        pz = Symbol('pz')

        coefs = pp.blade_coefs([self.e_0, self.e_1, self.e_2, self.e_3])
        p = coefs[0] + px * coefs[1] + py * coefs[2] + pz * coefs[3]
        self.assertEquals(p.subs({px: x1, py: y1, pz: z1}), S.Zero)
        self.assertEquals(p.subs({px: x2, py: y2, pz: z2}), S.Zero)
        self.assertEquals(p.subs({px: x3, py: y3, pz: z3}), S.Zero)

        coefs = pm.blade_coefs([self.e_0, self.e_1, self.e_2, self.e_3])
        p = coefs[0] + px * coefs[1] + py * coefs[2] + pz * coefs[3]
        self.assertEquals(p.subs({px: x1, py: y1, pz: z1}), S.Zero)
        self.assertEquals(p.subs({px: x2, py: y2, pz: z2}), S.Zero)
        self.assertEquals(p.subs({px: x3, py: y3, pz: z3}), S.Zero)

    def test_geometry_incidence_join_and_plane_side(self):
        """
        Join and plane side.
        """
        P1 = self.point(1, 0, 0)
        P2 = self.point(0, 1, 0)
        P3 = self.point(0, 0, 1)
        pp = Jinv(J(P1) ^ J(P2) ^ J(P3))

        px = Symbol('px')
        py = Symbol('py')
        pz = Symbol('pz')
        coefs = pp.blade_coefs([self.e_0, self.e_1, self.e_2, self.e_3])
        p = coefs[0] + px * coefs[1] + py * coefs[2] + pz * coefs[3]

        self.assertTrue(p.subs({px: Rational(1, 3) - 0.01, py: Rational(1, 3) - 0.01, pz: Rational(1, 3) - 0.01}) < 0.0)
        self.assertTrue(p.subs({px: Rational(1, 3) + 0.01, py: Rational(1, 3) + 0.01, pz: Rational(1, 3) + 0.01}) > 0.0)
