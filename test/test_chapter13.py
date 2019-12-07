from .test_utils import TestCase

from sympy import Symbol, S, sqrt
from galgebra.ga import Ga

USE_DIAGONAL_G = True


class TestChapter13(TestCase):

    def setUp(self):
        """
        Initialize CGA.
        """
        if USE_DIAGONAL_G:
            # This is way faster with galgebra but o and inf can't be printed...
            g = [
                [ 1,  0,  0,  0,  0],
                [ 0,  1,  0,  0,  0],
                [ 0,  0,  1,  0,  0],
                [ 0,  0,  0,  1,  0],
                [ 0,  0,  0,  0, -1],
            ]

            ga, e, e_1, e_2, e_3, eb = Ga.build('e e1 e2 e3 eb', g=g)
            o = S.Half * (e + eb)
            inf = eb - e

        else:
            g = [
                [ 0,  0,  0,  0, -1],
                [ 0,  1,  0,  0,  0],
                [ 0,  0,  1,  0,  0],
                [ 0,  0,  0,  1,  0],
                [-1,  0,  0,  0,  0],
            ]

            ga, o, e_1, e_2, e_3, inf = Ga.build('o e1 e2 e3 inf', g=g)

        self.ga = ga
        self.o = o
        self.e_1 = e_1
        self.e_2 = e_2
        self.e_3 = e_3
        self.inf = inf

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

    def symbol_vector(self, name='n', normalize=False):
        nx = Symbol(name + 'x', real=True)
        ny = Symbol(name + 'y', real=True)
        nz = Symbol(name + 'z', real=True)
        return self.vector(nx, ny, nz) / (sqrt(nx * nx + ny * ny + nz * nz) if normalize else S.One)

    def test_13_2_2(self):
        """
        Proper Euclidean motions as even versors : Translations.
        """
        delta_1 = Symbol('delta_1', real=True)
        delta_2 = Symbol('delta_2', real=True)
        n = self.symbol_vector('n', normalize=True)

        Tt = self.dual_plane(n, delta_2) * self.dual_plane(n, delta_1)

        # Even versor
        self.assertEqual(Tt * Tt.rev(), 1)

        # Translation
        r = Tt * self.o * Tt.inv()
        t = self.point(1, 2 * (delta_2 - delta_1) * n)
        self.assertEqual(r, t)

        # TODO: This exponential isn't available in galgebra
        #Te = (-t * self.inf * S.Half).exp()
        #self.assertEqual(Te * Te.rev(), 1)
        #self.assertEqual(Te, Tt)

    def test_13_2_3(self):
        """
        Proper Euclidean motions as even versors : Rotations in the origin.
        """
        n0 = self.symbol_vector('n0')
        n1 = self.symbol_vector('n1')
        R = n1 * n0     # 2 reflections
        Rinv = R.inv()

        p = self.symbol_vector('p')

        p0 = R * self.point(1, p) * Rinv
        p1 = self.point(1, R * p * Rinv)
        self.assertEqual(p0, p1)
