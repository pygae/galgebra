import unittest
from functools import reduce
from sympy import simplify, sqrt, Rational, Symbol
from galgebra.ga import Ga
from galgebra.mv import Mv
from galgebra.printer import Format, xtex
Format()
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

    def assertNotEquals(self, first, second, msg=None):
        """
        Compare two expressions are equals.
        """

        if isinstance(first, Mv):
            first = first.obj

        if isinstance(second, Mv):
            second = second.obj

        diff = simplify(first - second)

        self.assertTrue(diff != 0, "\n%s\n!=\n%s\n%s" % (first, second, diff))

    def test11_4(self):
        """
        All planes are 3-blades.
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        p0 = Symbol('p0')
        q0 = Symbol('q0')
        r0 = Symbol('r0')

        p = GA.mv((p0, Symbol('p1'), Symbol('p2'), Symbol('p3')), 'vector')
        q = GA.mv((q0, Symbol('q1'), Symbol('q2'), Symbol('q3')), 'vector')
        r = GA.mv((r0, Symbol('r1'), Symbol('r2'), Symbol('r3')), 'vector')

        p_inf = p.subs(p0,0)
        q_inf = q.subs(q0,0)
        r_inf = r.subs(r0,0)

        p = p.subs(p0,1)
        q = q.subs(q0,1)
        r = r.subs(r0,1)

        self.assertEquals(p ^ q ^ r, p ^ (q - p) ^ (r - p))
        self.assertEquals(p ^ q ^ r, p ^ (q_inf - p_inf) ^ (r_inf - p_inf))
        self.assertEquals(p ^ q ^ r, ((p + q + r) / 3) ^ ((p ^ q) + (q ^ r) + (r ^ p)))

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
            Ip_inv = Ip.inv()
            Ir = e_0_inv < Ip
            Ir_inv = Ir.inv()
            self.assertEquals(Ip, e_0 ^ Ir)
            self.assertEquals(Ip, e_0 * Ir)

            p = GA.mv([1] + [Symbol('p%d' % i) for i in range(1, GA.n)], 'vector')

            v = [
                GA.mv([0] + [Symbol('q%d' % i) for i in range(1, GA.n)], 'vector'),
                GA.mv([0] + [Symbol('r%d' % i) for i in range(1, GA.n)], 'vector'),
                GA.mv([0] + [Symbol('s%d' % i) for i in range(1, GA.n)], 'vector'),
                GA.mv([0] + [Symbol('t%d' % i) for i in range(1, GA.n)], 'vector'),
            ]

            # We test available finite k-flats
            for k in range(1, GA.n):
                A = reduce(Mv.__xor__, v[:k])
                X = (p ^ A)
                self.assertNotEquals(X, 0)
                M = e_0_inv < (e_0 ^ X)
                # Very slow
                d = (e_0_inv < (e_0 ^ X)) / (e_0_inv < X)
                print('d =',d.Fmt(3))
                #d_inv = d.inv()

                def hat(A):
                    return ((-1) ** A.pure_grade()) * A

                self.assertEquals(hat(A < Ir_inv), ((-1) ** (GA.n - 1)) * (hat(A) < Ir_inv))

                Xd = (p ^ A).dual()
                self.assertEquals(Xd, (p ^ A) < Ip_inv)
                self.assertEquals(Xd, p < (A < Ip_inv))
                self.assertEquals(Xd, p < ((A < Ir_inv) * e_0_inv))
                self.assertEquals(Xd, hat(A < Ir_inv) - e_0_inv * (p < hat(A < Ir_inv)))
                # Very slow
                #self.assertEquals(Xd, hat(A < Ir_inv) + e_0_inv * hat(M < Ir_inv))
                #self.assertEquals(Xd, (e_0_inv - d_inv) * hat(M < Ir_inv))
        xtex()

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
        self.assertEquals(p ^ q, p ^ d)
        self.assertEquals(p ^ q, q ^ d)

        r = e_0 + 2 * (e_2 - e_1)
        e = 2 * (e_2 - e_1)
        s = e_0 + 3 * (e_2 - e_1)
        t = 2 * (e_0 + e_2)
        self.assertEquals(2 * (p ^ q), p ^ e)
        self.assertEquals(2 * (p ^ q), p ^ t)

    def test11_12_2_1(self):
        """
        Let an orthonormal coordinate system be given in 3-dimensional Euclidean space.
        Compute the support vector of the line with direction u = e1 + 2 e2 - e3, through the point p = e1 - 3 e2.
        What is the distance of the line to the origin ?
        """
        GA, e_0, e_1, e_2, e_3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        L = (e_1 + 2 * e_2 - e_3) ^ (e_0 + e_1 - 3 * e_2)
        d = (e_0.inv() < (e_0 ^ L)) / (e_0.inv() < L)
        self.assertEquals(d.norm(), sqrt(Rational(35, 6)))

    """
    def test11_12_2_2(self):

        #Convert the line of the previous exercise into a parametric equation x = p + t * u; express t as a function of x
        #for a point x on the line.

        GA, e0, e1, e2, e3 = Ga.build("e*0|1|2|3", g='-1 0 0 0, 0 1 0 0, 0 0 1 0, 0 0 0 1')

        u = e1 + 2*e2 - e3
        p = e0 + e1 - 3*e2
        L = u ^ p

        t = Symbol('t')
        x1 = Symbol('x_1')
        x2 = Symbol('x_2')
        x3 = Symbol('x_3')
        x = e0 + x1 * e1 + x2 * e2 + x3 * e3

        # x(t)
        x_t = (e0 + e1 - 3 * e2) + t * (e1 + 2 * e2 - e3)

        d = e1 + 2 * e2 - e3
        d_inv = d.inv()

        # t(x)
        #t_x = (x - (e_0 + e_1 - 3 * e_2)) / (e_1 + 2 * e_2 - e_3)
        t_x = (x - (e0 + e1 - 3 * e2)) * d_inv

        for i in range(11):
            t_value = Rational(i, 10)
            x_value = x_t.subs(t,t_value)
            x_c = x_value.blade_coefs([e1, e2, e3])
            print(t_x.subs([x1,x2,x3],x_value.blade_coefs([e1, e2, e3])),' = ',t_value)
            self.assertEquals(x_value ^ L, 0)
            self.assertEquals(t_x.subs([x1, x2, x3], x_value.blade_coefs([e1, e2, e3])), t_value)
    """

if __name__ == '__main__':
    unittest.main()

