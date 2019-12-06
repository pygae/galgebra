from .test_utils import TestCase

from sympy import Symbol, S
from galgebra.ga import Ga


class TestChapter13(TestCase):

    def setUp(self):
        """
        Initialize CGA.
        """
        g = [
            [+0, +0, +0, +0, -1],
            [+0, +1, +0, +0, +0],
            [+0, +0, +1, +0, +0],
            [+0, +0, +0, +1, +0],
            [-1, +0, +0, +1, +0],
        ]

        self.ga, self.o, self.e_1, self.e_2, self.e_3, self.inf = Ga.build('o e1 e2 e3 inf', g=g)

    def vector(self, x, y, z):
        """
        Make a vector.
        """
        return x * self.e_1 + y * self.e_2 + z * self.e_3

    def point(self, alpha, v):
        """
        Make a point.
        """
        return alpha * (self.o + v + S.Half * v * v * self.inf)

    def dual_plane(self, n, delta):
        """
        Make a dual plane.
        """
        return n + delta * self.inf

    def dual_sphere(self, alpha, c, r):
        """
        Make a dual sphere.
        """
        return alpha * (c - S.Half * r * r * self.inf)

    def dual_im_sphere(self, alpha, c, r):
        """
        Make a dual imaginary sphere.
        """
        return alpha * (c + S.Half * r * r * self.inf)

    def test_13_1_2(self):
        """
        Points as null vectors.
        """
        p = self.point(1, self.vector(Symbol('px', real=True), Symbol('py', real=True), Symbol('pz', real=True)))
        q = self.point(1, self.vector(Symbol('qx', real=True), Symbol('qy', real=True), Symbol('qz', real=True)))

        self.assertEqual(p | q, -S.Half * q * q + p | q - S.Half * p * p)
        self.assertEqual(p | q, -S.Half * (q - p) * (q - p))

    def test_13_1_3(self):
        """
        General vectors represent dual planes and spheres.
        """
        alpha = Symbol('alpha', real=True)
        p = self.point(alpha, self.vector(Symbol('px', real=True), Symbol('py', real=True), Symbol('pz', real=True)))
        self.assertEqual(p | p, S.Zero)
        self.assertEqual(self.inf | p, -alpha)

        # Dual plane
        nx = Symbol('nx', real=True)
        ny = Symbol('ny', real=True)
        nz = Symbol('nz', real=True)
        p = self.dual_plane(self.vector(nx, ny, nz), Symbol('delta', real=True))
        self.assertEqual(p | p, nx * nx + ny * ny + nz * nz)
        self.assertEqual(self.inf | p, S.Zero)

        # Dual sphere
        cx = Symbol('cx', real=True)
        cy = Symbol('cy', real=True)
        cz = Symbol('cz', real=True)
        r = Symbol('r', real=True)

        c = self.point(1, self.vector(cx, cy, cz))
        self.assertEqual(c * c, S.Zero)
        self.assertEqual(-self.inf | c, S.One)

        s = self.dual_sphere(alpha, c, r)
        self.assertEqual(s | s, alpha * alpha * r * r)
        self.assertEqual(-self.inf | s, alpha)

        im_s = self.dual_im_sphere(alpha, c, r)
        self.assertEqual(im_s | im_s, -alpha * alpha * r * r)
        self.assertEqual(-self.inf | im_s, alpha)
