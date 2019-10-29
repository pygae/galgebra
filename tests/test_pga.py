import unittest

from sympy import acos, asin, pi, Rational, S, simplify, sqrt, Symbol, symbols

from ga import Ga
from mv import Mv, J, Jinv, com


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

    def norm(self, X):
        assert len(self.PGA.grade_decomposition(X)) == 1
        squared_norm = X | X
        if not squared_norm.is_scalar():
            raise ValueError("X | X isn't a scalar")
        squared_norm = squared_norm.scalar()
        if squared_norm == S.Zero:
            raise ValueError("X k-vector is null")
        return sqrt(abs(squared_norm))

    def norm2(self, X):
        assert len(self.PGA.grade_decomposition(X)) == 1
        squared_norm = X | X
        if not squared_norm.is_scalar():
            raise ValueError("X | X isn't a scalar")
        squared_norm = squared_norm.scalar()
        return squared_norm

    def ideal_norm(self, X):
        assert len(self.PGA.grade_decomposition(X)) == 1
        squared_norm = 0
        for c in X.blade_coefs():
            squared_norm += c * c
        if squared_norm == S.Zero:
            raise ValueError("X k-vector is null")
        return sqrt(squared_norm)

    def test_J(self):
        """
        Check we can join and meet using J and Jinv.
        """
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

        P = p1 ^ p2 ^ p3
        self.assertProjEquals(P, self.point(x, y, z))
        self.assertProjEquals(P, p1 ^ p3 ^ p2)
        self.assertProjEquals(P, p2 ^ p1 ^ p3)
        self.assertProjEquals(P, p2 ^ p3 ^ p1)
        self.assertProjEquals(P, p3 ^ p1 ^ p2)
        self.assertProjEquals(P, p3 ^ p2 ^ p1)

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
        #self.assertProjEquals(p123 ^ p124 ^ p314, P1)
        #self.assertProjEquals(p123 ^ p124 ^ p234, P2)
        #self.assertProjEquals(p123 ^ p234 ^ p314, P3)
        #self.assertProjEquals(p124 ^ p234 ^ p314, P4)

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
        Points join into lines.
        """
        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        x2, y2, z2 = symbols('x2 y2 z2')
        P2 = self.point(x2, y2, z2)

        lp = Jinv(J(P1) ^ J(P2))
        lm = Jinv(J(P2) ^ J(P1))
        self.assertEquals(lp, -lm)

        # TODO: add more tests...

    def test_metric_distance_of_points(self):
        """
        We can measure distance between normalized points using the joining line norm.
        """
        x0, y0, z0 = symbols('x0 y0 z0')
        P0 = self.point(x0, y0, z0)

        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        d = sqrt(abs((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2))

        lp = Jinv(J(P0) ^ J(P1))
        self.assertEquals(self.norm(lp), d)

        lm = Jinv(J(P1) ^ J(P0))
        self.assertEquals(self.norm(lm), d)

    def test_metric_angle_of_intersecting_planes(self):
        """
        We can measure the angle between normalized planes using the meeting line norm.
        """
        x = Symbol('x')
        y = Symbol('y')
        p1 = self.plane(1, 0, 0, -x)
        p2 = self.plane(0, 1, 0, -y)

        self.assertEquals(asin(self.norm(p1 ^ p2)), pi / 2)

    def test_metric_distance_between_parallel_planes(self):
        """
        We can measure the distance between two parallel and normalized planes using the meeting line norm.
        """
        nx, ny, nz, x, y = symbols('nx ny nz x y')

        n_norm = sqrt(nx * nx + ny * ny + nz * nz)
        p1 = self.plane(nx / n_norm, ny / n_norm, nz / n_norm, -x)
        p2 = self.plane(nx / n_norm, ny / n_norm, nz / n_norm, -y)

        self.assertEquals(self.ideal_norm(p1 ^ p2), sqrt((x - y)**2))

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

        # TODO: We could assume d0 - x0 is not null for simplifying further
        self.assertEquals(self.ideal_norm(P0 ^ p0), sqrt((d0 - x0)**2))

    def test_metric_oriented_distance_between_point_and_line(self):
        """
        We can measure the distance between a normalized point and a normalized line.
        """
        x0, y0, z0 = symbols('x0 y0 z0')
        P0 = self.point(x0, y0, z0)

        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        x2, y2, z2 = symbols('x2 y2 z2')
        P2 = self.point(x2, y2, z2)

        l0 = Jinv(J(P1) ^ J(P2))

        d = self.norm(Jinv(J(P0) ^ J(l0)))

        self.assertEquals(d.subs({x0: 2, y0: 0, z0: 0, x1: 0, y1: 2, z1: 0, x2: 1, y2: 2, z2: 0}), 2)
        self.assertEquals(d.subs({x0: 3, y0: 0, z0: 0, x1: 0, y1: 2, z1: 0, x2: 1, y2: 2, z2: 0}), 2)
        self.assertEquals(d.subs({x0: 2, y0: 0, z0: 0, x1: 0, y1: 1, z1: 0, x2: 1, y2: 1, z2: 0}), 1)
        self.assertEquals(d.subs({x0: 2, y0: 0, z0: 0, x1: 0, y1: 2, z1: 0, x2: 2, y2: 0, z2: 0}), 0)

    def test_metric_common_normal_line(self):
        """
        We can find the common normal line of two normalized lines.
        """
        x0, y0, z0 = symbols('x0 y0 z0')
        P0 = self.point(x0, y0, z0)

        x1, y1, z1 = symbols('x1 y1 z1')
        P1 = self.point(x1, y1, z1)

        x2, y2, z2 = symbols('x2 y2 z2')
        P2 = self.point(x2, y2, z2)

        x3, y3, z3 = symbols('x3 y3 z3')
        P3 = self.point(x3, y3, z3)

        l0 = Jinv(J(P0) ^ J(P1))
        l1 = Jinv(J(P2) ^ J(P3))
        ln = com(l0, l1)

        ln /= self.norm(ln)
        l0 /= self.norm(l0)
        l1 /= self.norm(l1)

        self.assertEquals(acos(l0 | ln), S.Half * pi)
        self.assertEquals(acos(l1 | ln), S.Half * pi)

    def test_metric_angle_between_lines(self):
        """
        We can measure the angle between to normalized lines.
        """
        x1, y1 = symbols('x1 y1', real=True)
        P0 = self.point(0, 0, 0)
        P1 = self.point(x1, y1, 0)
        P2 = self.point(-y1, x1, 0)

        l0 = Jinv(J(P0) ^ J(P1))    # TODO: this feels weird... but ganja does the same
        l1 = Jinv(J(P2) ^ J(P0))    #

        l0 /= self.norm(l0)
        l1 /= self.norm(l1)

        self.assertEquals(acos(l0 | l1), pi / 2)
