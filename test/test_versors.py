import unittest
from galgebra.ga import Ga
from sympy import S


class TestVersors(unittest.TestCase):
    def setUp(self):
        # Set up different geometric algebras for testing

        # G(3,0) - Euclidean 3-space
        self.g3d = Ga("e1 e2 e3", g=[1, 1, 1])

        # G(1,3) - Spacetime algebra
        self.sta = Ga("e0 e1 e2 e3", g=[1, -1, -1, -1])

        # Initialize basis vectors
        e1, e2, e3 = self.g3d.mv()
        e0, e1, e2, e3 = self.sta.mv()

    def test_basis_vectors_are_versors(self):
        """Individual basis vectors should be versors"""
        e1, e2, e3 = self.g3d.mv()
        self.assertTrue(e1.is_versor())
        self.assertTrue(e2.is_versor())
        self.assertTrue(e3.is_versor())

    def test_mixed_grades_not_versor(self):
        """A sum of different grades cannot be a versor"""
        e1, e2, e3 = self.g3d.mv()
        mixed = e1 + e2 * e3  # vector + bivector
        self.assertFalse(mixed.is_versor())

    def test_dorst_counterexample(self):
        """Test Greg1950's counterexample from the spacetime algebra"""
        e0, e1, e2, e3 = self.sta.mv()
        V = S.One + self.sta.I()  # 1 + I
        self.assertFalse(V.is_versor())

    def test_rotors_are_versors(self):
        """Test that rotors (even-grade versors) are properly detected"""
        e1, e2, e3 = self.g3d.mv()
        rotor = S.One + e1 * e2  # 1 + e1e2
        self.assertTrue(rotor.is_versor())

    def test_null_is_not_versor(self):
        """Test that null multivectors are not versors"""
        e1, e2, e3 = self.g3d.mv()
        null_mv = e1 + e1 * self.g3d.I()  # v + vI (squares to 0)
        self.assertFalse(null_mv.is_versor())
