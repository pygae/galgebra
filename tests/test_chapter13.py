from .test_utils import TestCase

from sympy import Symbol, S
from galgebra.ga import Ga


class TestChapter13(TestCase):

    def test_13_1_2(self):
        """
        Points as null vectors.
        """
        g = [
            [+0, +0, +0, +0, -1],
            [+0, +1, +0, +0, +0],
            [+0, +0, +1, +0, +0],
            [+0, +0, +0, +1, +0],
            [-1, +0, +0, +1, +0],
        ]

        GA, o, e_1, e_2, e_3, inf = Ga.build('o e1 e2 e3 inf', g=g)

        def vector(x, y, z):
            return x * e_1 + y * e_2 + z * e_3

        def point(v):
            return o + v + S.Half * v * v * inf

        p = point(vector(Symbol('px', real=True), Symbol('py', real=True), Symbol('pz', real=True)))
        q = point(vector(Symbol('qx', real=True), Symbol('qy', real=True), Symbol('qz', real=True)))

        self.assertEquals(p | q, -S.Half * q * q + p | q - S.Half * p * p)
        self.assertEquals(p | q, -S.Half * (q - p) * (q - p))

    def test_13_1_3(self):
        """
        General vectors represent dual planes and spheres.
        """
        g = [
            [+0, +0, +0, +0, -1],
            [+0, +1, +0, +0, +0],
            [+0, +0, +1, +0, +0],
            [+0, +0, +0, +1, +0],
            [-1, +0, +0, +1, +0],
        ]

        GA, o, e_1, e_2, e_3, inf = Ga.build('o e1 e2 e3 inf', g=g)

        def vector(x, y, z):
            return x * e_1 + y * e_2 + z * e_3

        # Point
        def point(alpha, v):
            return alpha * (o + v + S.Half * v * v * inf)

        alpha = Symbol('alpha', real=True)
        p = point(alpha, vector(Symbol('px', real=True), Symbol('py', real=True), Symbol('pz', real=True)))
        self.assertEquals(p | p, S.Zero)
        self.assertEquals(inf | p, -alpha)

        # Dual plane
        def dual_plane(n, delta):
            return n + delta * inf

        nx = Symbol('nx', real=True)
        ny = Symbol('ny', real=True)
        nz = Symbol('nz', real=True)
        p = dual_plane(vector(nx, ny, nz), Symbol('delta', real=True))
        self.assertEquals(p | p, nx * nx + ny * ny + nz * nz)
        self.assertEquals(inf | p, S.Zero)

        # Dual sphere
        def dual_sphere(alpha, c, r):
            return alpha * (c - S.Half * r * r * inf)

        def dual_im_sphere(alpha, c, r):
            return alpha * (c + S.Half * r * r * inf)

        cx = Symbol('cx', real=True)
        cy = Symbol('cy', real=True)
        cz = Symbol('cz', real=True)
        r = Symbol('r', real=True)

        c = point(1, vector(cx, cy, cz))
        self.assertEquals(c * c, S.Zero)
        self.assertEquals(-inf | c, S.One)

        s = dual_sphere(alpha, c, r)
        self.assertEquals(s | s, alpha * alpha * r * r)
        self.assertEquals(-inf | s, alpha)

        im_s = dual_im_sphere(alpha, c, r)
        self.assertEquals(im_s | im_s, -alpha * alpha * r * r)
        self.assertEquals(-inf | im_s, alpha)
