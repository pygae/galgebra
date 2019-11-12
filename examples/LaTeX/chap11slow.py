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
        """
        -1 0 0 0 0
         0 1 0 0 0
         0 0 1 0 0
         0 0 0 1 0
         0 0 0 0 1
        """

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

        Ga.dual_mode()

if __name__ == '__main__':
    unittest.main()

