from .test_utils import TestCase

from functools import reduce
from sympy import sqrt, Rational, Symbol
from galgebra.ga import Ga
from galgebra.mv import Mv


class TestChapter11(TestCase):

    def test11_4(self):
        """
        All planes are 3-blades.
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        p0 = Symbol('p0', real=True)
        q0 = Symbol('q0', real=True)
        r0 = Symbol('r0', real=True)

        p = GA.mv((p0, Symbol('p1', real=True), Symbol('p2', real=True), Symbol('p3', real=True)), 'vector')
        q = GA.mv((q0, Symbol('q1', real=True), Symbol('q2', real=True), Symbol('q3', real=True)), 'vector')
        r = GA.mv((r0, Symbol('r1', real=True), Symbol('r2', real=True), Symbol('r3', real=True)), 'vector')

        p_inf = p.subs({p0: 0})
        q_inf = q.subs({q0: 0})
        r_inf = r.subs({r0: 0})

        p = p.subs({p0: 1})
        q = q.subs({q0: 1})
        r = r.subs({r0: 1})

        self.assertEqual(p ^ q ^ r, p ^ (q - p) ^ (r - p))
        self.assertEqual(p ^ q ^ r, p ^ (q_inf - p_inf) ^ (r_inf - p_inf))
        self.assertEqual(p ^ q ^ r, ((p + q + r) / 3) ^ ((p ^ q) + (q ^ r) + (r ^ p)))

    def test11_6(self):
        """
        Dual representation.
        """
        Ga.dual_mode('Iinv+')

        GA_list = [
            Ga("e*0|1|2", g='-1 0 0, 0 1 0, 0 0 1'),
            Ga("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1'),
            Ga("e*0|1|2|3|5", g='-1 0 0 0 0, 0 1 0 0 0, 0 0 1 0 0, 0 0 0 1 0, 0 0 0 0 1'),
        ]

        for GA in GA_list:
            e_0 = GA.mv_basis[0]
            e_0_inv = e_0.inv()

            Ip = GA.I()
            Ip_inv = GA.I_inv()
            Ir = e_0_inv < Ip
            Ir_inv = Ir.inv()
            self.assertEqual(Ip, e_0 ^ Ir)
            self.assertEqual(Ip, e_0 * Ir)

            p = GA.mv([1] + [Symbol('p%d' % i, real=True) for i in range(1, GA.n)], 'vector')

            v = [
                GA.mv([0] + [Symbol('q%d' % i, real=True) for i in range(1, GA.n)], 'vector'),
                GA.mv([0] + [Symbol('r%d' % i, real=True) for i in range(1, GA.n)], 'vector'),
                GA.mv([0] + [Symbol('s%d' % i, real=True) for i in range(1, GA.n)], 'vector'),
                GA.mv([0] + [Symbol('t%d' % i, real=True) for i in range(1, GA.n)], 'vector'),
            ]

            # We test available finite k-flats
            for k in range(1, GA.n):
                A = reduce(Mv.__xor__, v[:k])
                X = (p ^ A)
                self.assertNotEqual(X, 0)
                M = e_0_inv < (e_0 ^ X)

                # Very slow
                # d = (e_0_inv < (e_0 ^ X)) / (e_0_inv < X)
                # d_inv = d.inv()

                def hat(A):
                    return ((-1) ** A.pure_grade()) * A

                self.assertEqual(hat(A < Ir_inv), ((-1) ** (GA.n - 1)) * (hat(A) < Ir_inv))

                Xd = (p ^ A).dual()
                self.assertEqual(Xd, (p ^ A) < Ip_inv)
                self.assertEqual(Xd, p < (A < Ip_inv))
                self.assertEqual(Xd, p < ((A < Ir_inv) * e_0_inv))
                self.assertEqual(Xd, hat(A < Ir_inv) - e_0_inv * (p < hat(A < Ir_inv)))
                # Very slow
                # self.assertEquals(Xd, hat(A < Ir_inv) + e_0_inv * hat(M < Ir_inv))
                # self.assertEquals(Xd, (e_0_inv - d_inv) * hat(M < Ir_inv))

        Ga.dual_mode()

    def test11_12_1(self):
        """
        Compute the 2-blades corresponding to the lines gives by the data below. Which of
        the lines are the same, considered as weighted oriented elements of geometry, which
        are the same as offset subspaces ?
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        p = e_0 + e_1
        q = e_0 + e_2
        d = e_2 - e_1
        self.assertEqual(p ^ q, p ^ d)
        self.assertEqual(p ^ q, q ^ d)

        r = e_0 + 2 * (e_2 - e_1)
        e = 2 * (e_2 - e_1)
        s = e_0 + 3 * (e_2 - e_1)
        t = 2 * (e_0 + e_2)
        self.assertEqual(2 * (p ^ q), p ^ e)
        self.assertEqual(2 * (p ^ q), p ^ t)

    def test11_12_2_1(self):
        """
        Let an orthonormal coordinate system be given in 3-dimensional Euclidean space.
        Compute the support vector of the line with direction u = e1 + 2 e2 - e3, through the point p = e1 - 3 e2.
        What is the distance of the line to the origin ?
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        L = (e_1 + 2 * e_2 - e_3) ^ (e_0 + e_1 - 3 * e_2)
        d = (e_0.inv() < (e_0 ^ L)) / (e_0.inv() < L)
        self.assertEqual(d.norm(), sqrt(Rational(35, 6)))

    def test11_12_2_2(self):
        """
        Convert the line of the previous exercise into a parametric equation x = p + t * u; express t as a function of x
        for a point x on the line.
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        u = e_1 + 2 * e_2 - e_3
        p = e_0 + e_1 - 3 * e_2
        L = u ^ p

        t = Symbol('t', real=True)
        x_1 = Symbol('x_1', real=True)
        x_2 = Symbol('x_2', real=True)
        x_3 = Symbol('x_3', real=True)
        x = e_0 + x_1 * e_1 + x_2 * e_2 + x_3 * e_3

        # x(t)
        x_t = (e_0 + e_1 - 3 * e_2) + t * (e_1 + 2 * e_2 - e_3)

        # t(x)
        t_x = (x - (e_0 + e_1 - 3 * e_2)) / (e_1 + 2 * e_2 - e_3)

        for i in range(11):
            t_value = Rational(i, 10)
            x_value = x_t.subs({t: t_value})
            self.assertEqual(x_value ^ L, 0)
            self.assertEqual(t_x.subs(list(zip([x_1, x_2, x_3], x_value.blade_coefs([e_1, e_2, e_3])))), t_value)
