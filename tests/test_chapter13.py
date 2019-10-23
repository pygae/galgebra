import unittest

from sympy import simplify, Symbol, S

from ga import Ga
from mv import Mv


class TestChapter13(unittest.TestCase):

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

        p = point(vector(Symbol('px'), Symbol('py'), Symbol('pz')))
        q = point(vector(Symbol('qx'), Symbol('qy'), Symbol('qz')))

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

        alpha = Symbol('alpha')
        p = point(alpha, vector(Symbol('px'), Symbol('py'), Symbol('pz')))
        self.assertEquals(p | p, S.Zero)
        self.assertEquals(inf | p, -alpha)

        # Dual plane
        def dual_plane(n, delta):
            return n + delta * inf

        nx = Symbol('nx')
        ny = Symbol('ny')
        nz = Symbol('nz')
        p = dual_plane(vector(nx, ny, nz), Symbol('delta'))
        self.assertEquals(p | p, nx * nx + ny * ny + nz * nz)
        self.assertEquals(inf | p, S.Zero)

        # Dual sphere
        def dual_sphere(alpha, c, r):
            return alpha * (c - S.Half * r * r * inf)

        def dual_im_sphere(alpha, c, r):
            return alpha * (c + S.Half * r * r * inf)

        cx = Symbol('cx')
        cy = Symbol('cy')
        cz = Symbol('cz')
        r = Symbol('r')

        c = point(1, vector(cx, cy, cz))
        self.assertEquals(c * c, S.Zero)
        self.assertEquals(-inf | c, S.One)

        s = dual_sphere(alpha, c, r)
        self.assertEquals(s | s, alpha * alpha * r * r)
        self.assertEquals(-inf | s, alpha)

        im_s = dual_im_sphere(alpha, c, r)
        self.assertEquals(im_s | im_s, -alpha * alpha * r * r)
        self.assertEquals(-inf | im_s, alpha)
