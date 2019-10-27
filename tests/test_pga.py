import unittest

from sympy import Abs, asin, pi, Rational, S, simplify, sqrt, Symbol, symbols

from ga import Ga
from mv import Mv, J, Jinv


class TestPGA(unittest.TestCase):

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
        self.e_01 = e_0 ^ e_1
        self.e_02 = e_0 ^ e_2
        self.e_03 = e_0 ^ e_3
        self.e_12 = e_1 ^ e_2
        self.e_31 = e_3 ^ e_1
        self.e_23 = e_2 ^ e_3
        self.e_032 = e_0 ^ e_3 ^ e_2
        self.e_013 = e_0 ^ e_1 ^ e_3
        self.e_021 = e_0 ^ e_2 ^ e_1
        self.e_123 = e_1 ^ e_2 ^ e_3
        self.e_0123 = e_0 ^ e_1 ^ e_2 ^ e_3

        self.e_13 = e_1 ^ e_3

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

    def line_norm(self, l, hint=None):
        assert hint == 'euclidean' or hint == 'ideal'
        if hint == 'euclidean':
            d, e, f = l.blade_coefs([self.e_12, self.e_13, self.e_23])
            return sqrt(d * d + e * e + f * f)
        elif hint == 'ideal':
            a, b, c = l.blade_coefs([self.e_01, self.e_02, self.e_03])
            return sqrt(a * a + b * b + c * c)

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

        print pp
        print pm

        a = Symbol('a')
        b = Symbol('b')
        c = Symbol('c')

        coefs = pp.blade_coefs([self.e_0, self.e_1, self.e_2, self.e_3])
        p = coefs[0] + a * coefs[1] + b * coefs[2] + c * coefs[3]
        self.assertEquals(p.subs({a: x1, b: y1, c: z1}), S.Zero)
        self.assertEquals(p.subs({a: x2, b: y2, c: z2}), S.Zero)
        self.assertEquals(p.subs({a: x3, b: y3, c: z3}), S.Zero)

        coefs = pm.blade_coefs([self.e_0, self.e_1, self.e_2, self.e_3])
        p = coefs[0] + a * coefs[1] + b * coefs[2] + c * coefs[3]
        self.assertEquals(p.subs({a: x1, b: y1, c: z1}), S.Zero)
        self.assertEquals(p.subs({a: x2, b: y2, c: z2}), S.Zero)
        self.assertEquals(p.subs({a: x3, b: y3, c: z3}), S.Zero)

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

    def test_geometry_incidence_points_join_into_lines(self):
        """
        Points join into lines, planes meet into lines.
        """
        PGA = self.PGA

        x0, y0, z0 = symbols('x0 y0 z0')
        P0 = self.point(x0, y0, z0)

        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        x2, y2, z2 = symbols('x2 y2 z2')
        P2 = self.point(x2, y2, z2)

        x3, y3, z3 = symbols('x3 y3 z3')
        P3 = self.point(x3, y3, z3)

        lp = Jinv(J(P1) ^ J(P2))
        lm = Jinv(J(P2) ^ J(P1))
        self.assertEquals(lp, -lm)

        p012 = Jinv(J(P0) ^ J(P1) ^ J(P2))
        p123 = Jinv(J(P1) ^ J(P2) ^ J(P3))
        l = p012 ^ p123

        print l / self.line_norm(l)

    def test_metric_distance_of_points(self):
        """
        We can measure distance between normalized points using the joining line norm.
        """
        x0, y0, z0 = symbols('x0 y0 z0')
        P0 = self.point(x0, y0, z0)

        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        d = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)

        lp = Jinv(J(P0) ^ J(P1))
        self.assertEquals(self.line_norm(lp, 'euclidean'), d)

        lm = Jinv(J(P1) ^ J(P0))
        self.assertEquals(self.line_norm(lm, 'euclidean'), d)

    def test_metric_angle_of_intersecting_planes(self):
        """
        We can measure the angle between normalized planes using the meeting line norm.
        """
        x = Symbol('x')
        y = Symbol('y')
        p1 = self.plane(1, 0, 0, -x)
        p2 = self.plane(0, 1, 0, -y)

        self.assertEquals(asin(self.line_norm(p1 ^ p2, 'euclidean')), pi / 2)

    def test_metric_distance_between_parallel_planes(self):
        """
        We can measure the distance between two parallel and normalized planes using the meeting line norm.
        """
        nx, ny, nz, x, y = symbols('nx ny nz x y')

        n_norm = sqrt(nx * nx + ny * ny + nz * nz)
        p1 = self.plane(nx / n_norm, ny / n_norm, nz / n_norm, -x)
        p2 = self.plane(nx / n_norm, ny / n_norm, nz / n_norm, -y)

        self.assertEquals(self.line_norm(p1 ^ p2, 'ideal'), sqrt((x - y)**2))

    def test_metric_oriented_distance_between_point_and_plane(self):
        """
        We can measure the distance between a normalized point and a normalized plane.
        """

        x0, y0, z0 = symbols('x0 y0 z0')
        P0 = self.point(x0, y0, z0)

        d0 = Symbol('d0')
        p0 = self.plane(1, 0, 0, -d0)

        self.assertEquals(Jinv(J(P0) ^ J(p0)), d0 - x0)
        self.assertEquals(Jinv(J(p0) ^ J(P0)), x0 - d0)

        self.assertEquals(P0 ^ p0, (d0 - x0) * self.e_0123)     # TODO: how can we use inf norm here ?
        self.assertEquals(p0 ^ P0, (x0 - d0) * self.e_0123)     #
