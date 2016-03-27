import unittest

from itertools import product, izip, count
from sympy import Symbol
from ga import Ga

class TestChapter3(unittest.TestCase):
    
    def test_3_2_2(self):
        """
        Computing the contraction explicitly.
        """
        
        (_g3d, e_1, e_2, e_3) = Ga.build('e*1|2|3')
                
        A_s = Symbol('A_s')
        A_1 = Symbol('A_1')
        A_2 = Symbol('A_2')
        A_3 = Symbol('A_3')
        A_12 = Symbol('A_12')
        A_13 = Symbol('A_13')
        A_23 = Symbol('A_23')
        A_123 = Symbol('A_123')
        B_s = Symbol('B_s')
        B_1 = Symbol('B_1')
        B_2 = Symbol('B_2')
        B_3 = Symbol('B_3')
        B_12 = Symbol('B_12')
        B_13 = Symbol('B_13')
        B_23 = Symbol('B_23')
        B_123 = Symbol('B_123')
        C_s = Symbol('C_s')
        C_1 = Symbol('C_1')
        C_2 = Symbol('C_2')
        C_3 = Symbol('C_3')
        C_12 = Symbol('C_12')
        C_13 = Symbol('C_13')
        C_23 = Symbol('C_23')
        C_123 = Symbol('C_123')        
        
        A_by_grades = [
             _g3d.mv(A_s),
             A_1 * e_1 + A_2 * e_2 + A_3 * e_3,
             A_12 * (e_1 ^ e_2) + A_13 * (e_1 ^ e_3) + A_23 * (e_2 ^ e_3),
             A_123 * (e_1 ^ e_2 ^ e_3),
        ]
        
        B_by_grades = [
             _g3d.mv(B_s),
             B_1 * e_1 + B_2 * e_2 + B_3 * e_3,
             B_12 * (e_1 ^ e_2) + B_13 * (e_1 ^ e_3) + B_23 * (e_2 ^ e_3),
             B_123 * (e_1 ^ e_2 ^ e_3),
        ]
        
        C_by_grades = [
             _g3d.mv(C_s),
             C_1 * e_1 + C_2 * e_2 + C_3 * e_3,
             C_12 * (e_1 ^ e_2) + C_13 * (e_1 ^ e_3) + C_23 * (e_2 ^ e_3),
             C_123 * (e_1 ^ e_2 ^ e_3),
        ]
        
        self.assertTrue(A_by_grades[1].is_vector())
        self.assertTrue(B_by_grades[1].is_vector())
        self.assertTrue(C_by_grades[1].is_vector())
        
        for grade, A in izip(count(), A_by_grades):
            self.assertTrue(A.is_blade() and A.pure_grade() == grade)
        
        for grade, B in izip(count(), B_by_grades):
            self.assertTrue(B.is_blade() and B.pure_grade() == grade)
        
        for grade, C in izip(count(), C_by_grades):
            self.assertTrue(C.is_blade() and C.pure_grade() == grade)
        
        # scalar and blades of various grades
        A = A_by_grades[0]
        for B in B_by_grades:
            self.assertTrue((A < B) == (A * B))
            
        A = A_by_grades[0]
        for B in B_by_grades:
            self.assertTrue((B < A) == 0 if B.pure_grade() > 0 else B_s)
        
        # vectors
        A = A_by_grades[1]
        B = B_by_grades[1]        
        self.assertTrue((A < B) == (A | B))
        
        # vector and the outer product of 2 blades of various grades (scalars, vectors, 2-vectors...)
        A = A_by_grades[1]
        for B, C in product(B_by_grades, C_by_grades):
            self.assertTrue((A < (B ^ C)) == (((A < B) ^ C) + (-1)**B.pure_grade() * (B ^ (A < C))))
            
        # vector and the outer product of 2 blades of various grades (scalars, vectors, 2-vectors...)
        for A, B, C in product(A_by_grades, B_by_grades, C_by_grades):
            self.assertTrue(((A ^ B) < C) == (A < (B < C)))
    

if __name__ == '__main__':

    unittest.main()
