import unittest
import itertools

from sympy import Symbol, Matrix, solve, solve_poly_system, cos, sin
from ga import Ga
from mv import Mv

class TestChapter2(unittest.TestCase):

    def test_2_12_1_1(self):
        """
        Compute the outer products of the following 3-space expressions,
        giving the result relative to the basis {1, e_1, e_2, e_3, e_1^e_2, e_1^e_3, e_2^e_3, e_1^e_2^e_3}.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')
        self.assertTrue((e_1 + e_2) ^ (e_1 + e_3) == (-e_1 ^ e_2) + (e_1 ^ e_3) + (e_2 ^ e_3))
        self.assertTrue((e_1 + e_2 + e_3) ^ (2*e_1) == -2*(e_1 ^ e_2) - 2*(e_1 ^ e_3))
        self.assertTrue((e_1 - e_2) ^ (e_1 - e_3) == (e_1 ^ e_2) - (e_1 ^ e_3) + (e_2 ^ e_3))
        self.assertTrue((e_1 + e_2) ^ (0.5*e_1 + 2*e_2 + 3*e_3) == 1.5*(e_1 ^ e_2) + 3*(e_1 ^ e_3) + 3*(e_2 ^ e_3))
        self.assertTrue((e_1 ^ e_2) ^ (e_1 + e_3) == (e_1 ^ e_2 ^ e_3))
        self.assertTrue((e_1 + e_2) ^ ((e_1 ^ e_2) + (e_2 ^ e_3)) == (e_1 ^ e_2 ^ e_3))


    def test2_12_1_2(self):
        """
        Given the 2-blade B = e_1 ^ (e_2 - e_3) that represents a plane,
        determine if each of the following vectors lies in that plane.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')
        B = e_1 ^ (e_2 - e_3)
        self.assertTrue(e_1 ^ B == 0)
        self.assertFalse((e_1 + e_2) ^ B == 0)
        self.assertFalse((e_1 + e_2 + e_3) ^ B == 0)
        self.assertTrue((2*e_1 - e_2 + e_3) ^ B == 0)


    def test2_12_1_3(self):
        """
        What is the area of the parallelogram spanned by the vectors a = e_1 + 2*e_2
        and b = -e_1 - e_2 (relative to the area of e_1 ^ e_2) ?
        """
        (_g3d, e_1, e_2, _e_3) = Ga.build('e*1|2|3')
        a = e_1 + 2*e_2
        b = -e_1 - e_2
        B = a ^ b
        self.assertTrue(B == 1 * (e_1 ^ e_2))


    def test2_12_1_4(self):
        """
        Compute the intersection of the non-homogeneous line L with position vector e_1
        and direction vector e_2, and the line M with position vector e_2 and direction vector
        (e_1 + e_2), using 2-blades.
        """
        (_g2d, e_1, e_2) = Ga.build('e*1|2')

        # x is defined in the basis {e_1, e_2}
        a = Symbol('a')
        b = Symbol('b')
        x = a * e_1 + b * e_2

        # x lies on L and M (eq. L == 0 and M == 0)
        L = (x ^ e_2) - (e_1 ^ e_2)
        M = (x ^ (e_1 + e_2)) - (e_2 ^ (e_1 + e_2))

        # Solve the linear system
        R = solve([L, M], a, b)

        # Replace symbols
        x = x.subs(R)

        self.assertTrue(x == e_1 + 2*e_2)


    def test_2_12_1_5(self):
        """
        Compute (2 + 3 * e_3) ^ (e_1 + (e_2 ^ e_3) using the grade based defining equations of section 2.9.4.
        """
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')

        self.assertTrue((2 + 3 * e_3) ^ (e_1 + (e_2 ^ e_3)) == (2 * e_1 + 2 * (e_2 ^ e_3) + 3 * (e_3 ^ e_1)))


    def test2_12_2_1(self):
        """
        In R2 with Euclidean metric, choose an orthonormal basis {e_1, e_2} in the plane of a and b such that e1 is parallel to a.
        Write x = a * e_1 and y = b * (cos(t) * e_1 + sin(t) * e_2), whete t is the angle from a to b.
        Evaluate the outer product. What is the geometrical interpretation ?
        """
        (_g2d, e_1, e_2) = Ga.build('e*1|2', g='1 0, 0 1')

        # TODO: use alpha, beta and theta instead of a, b and t (it crashes sympy)
        a = Symbol('a')
        b = Symbol('b')
        t = Symbol('t')
        x = a * e_1
        y = b * (cos(t) * e_1 + sin(t) * e_2)
        B = x ^ y
        self.assertTrue(B == (a * b * sin(t) * (e_1 ^ e_2)))

        # Retrieve the parallelogram area from the 2-vector
        area = B.norm()
        self.assertTrue(area == (a * b * sin(t)))

        # Compute the parallelogram area using the determinant
        x = [a, 0]
        y = [b * cos(t), b * sin(t)]
        area = Matrix([x, y]).det()
        self.assertTrue(area == (a * b * sin(t)))


    def test2_12_2_2(self):

        (_g2d, e_1, e_2) = Ga.build('e*1|2', g='1 0, 0 1')

        a = Symbol('a')
        b = Symbol('b')
        t = Symbol('t')
        x = a * e_1
        y = b * (cos(t) * e_1 + sin(t) * e_2)

        x_1 = x.blade_coefs([e_1])[0]
        x_2 = x.blade_coefs([e_2])[0]
        y_1 = y.blade_coefs([e_1])[0]
        y_2 = y.blade_coefs([e_2])[0]

        cross = x_1 * y_2 - y_1 * x_2

        self.assertTrue(cross == (a * b * sin(t)))


    def test2_12_2_4(self):
        """
        Solve the linear system of equation using the outer product.
        """
        (g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')

        a_1 = Symbol('a_1')
        a_2 = Symbol('a_2')
        a_3 = Symbol('a_3')
        b_1 = Symbol('b_1')
        b_2 = Symbol('b_2')
        b_3 = Symbol('b_3')
        c_1 = Symbol('c_1')
        c_2 = Symbol('c_2')
        c_3 = Symbol('c_3')
        a = a_1 * e_1 + a_2 * e_2 + a_3 * e_3
        b = b_1 * e_1 + b_2 * e_2 + b_3 * e_3
        c = c_1 * e_1 + c_2 * e_2 + c_3 * e_3

        x_1 = Symbol('x_1')
        x_2 = Symbol('x_2')
        x_3 = Symbol('x_3')
        x = x_1 * a + x_2 * b + x_3 * c

        # Solve x_1
        self.assertTrue((x ^ a) == (x_2 * (b ^ a) + x_3 * (c ^ a)))
        self.assertTrue((x ^ a ^ b) == x_3 * (c ^ a ^ b))
        self.assertTrue((x ^ a ^ b) * (c ^ a ^ b).inv() == Mv(x_3, 'scalar', ga=g3d))

        # Solve x_2
        self.assertTrue((x ^ b) == (x_1 * (a ^ b) + x_3 * (c ^ b)))
        self.assertTrue((x ^ b ^ c) == x_1 * (a ^ b ^ c))
        self.assertTrue((x ^ b ^ c) * (a ^ b ^ c).inv() == Mv(x_1, 'scalar', ga=g3d))

        # Solve x_3
        self.assertTrue((x ^ c) == (x_1 * (a ^ c) + x_2 * (b ^ c)))
        self.assertTrue((x ^ c ^ a) == x_2 * (b ^ c ^ a))
        self.assertTrue((x ^ c ^ a) * (b ^ c ^ a).inv() == Mv(x_2, 'scalar', ga=g3d))


    def test2_12_2_5(self):
        """
        Consider R4 with basis {e_1, e_2, e_3, e_4}. Show that the 2-vector B = (e_1 ^ e_2) + (e_3 ^ e_4)
        is not a 2-blade (i.e., it cannot be written as the outer product of two vectors).
        """
        (_g4d, e_1, e_2, e_3, e_4) = Ga.build('e*1|2|3|4')

        # B
        B = (e_1 ^ e_2) + (e_3 ^ e_4)

        # C is the product of a and b vectors
        a_1 = Symbol('a_1')
        a_2 = Symbol('a_2')
        a_3 = Symbol('a_3')
        a_4 = Symbol('a_4')
        a = a_1 * e_1 + a_2 * e_2 + a_3 * e_3 + a_4 * e_4

        b_1 = Symbol('b_1')
        b_2 = Symbol('b_2')
        b_3 = Symbol('b_3')
        b_4 = Symbol('b_4')
        b = b_1 * e_1 + b_2 * e_2 + b_3 * e_3 + b_4 * e_4

        C = a ^ b

        # other coefficients are null
        blades = [
            e_1 ^ e_2, e_1 ^ e_3, e_1 ^ e_4, e_2 ^ e_3, e_2 ^ e_4, e_3 ^ e_4,
        ]

        C_coefs = C.blade_coefs(blades)
        B_coefs = B.blade_coefs(blades)

        # try to solve the system and show there is no solution
        system = [
            (C_coef) - (B_coef) for C_coef, B_coef in zip(C_coefs, B_coefs)
        ]

        unknowns = [
            a_1, a_2, a_3, a_4, b_1, b_2, b_3, b_4
        ]

        # TODO: use solve if sympy fix it
        result = solve_poly_system(system, unknowns)
        self.assertTrue(result is None)


    def test2_12_2_9(self):
        """
        Prove Ak ^ Bl = (-1**kl) Bl ^ Ak.
        """

        # very slow if bigger
        dimension = 5

        # create the vector space
        # e*1|2|3|4|5 if dimension == 5
        # TODO: don't fix the vector space dimension
        ga = Ga('e*%s' % '|'.join(str(i) for i in range(1, dimension+1)))

        # an helper for building multivectors
        # Ak = a1 ^ a2 ^ ... ^ ak if _build('a', k)
        def _build(name, grade):
            M = ga.mv('%s1' % name, 'vector')
            for k in range(2, grade+1):
                M = M ^ ga.mv('%s%d' % (name, k), 'vector')
            return M

        # prove for k in [1, 5] and l in [1, 5]
        for k, l in itertools.product(range(1, dimension+1), range(1, dimension+1)):
            Ak = _build('a', k)
            Bl = _build('b', l)
            self.assertTrue((Ak ^ Bl) == (-1)**(k*l) * (Bl ^ Ak))


if __name__ == '__main__':

    unittest.main()

