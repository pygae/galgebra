"""
Validation of galgebra against Russell Goyder's geometric algebra cheat sheet.

Cheat sheet: https://russellgoyder.github.io/geometric-algebra-cheat-sheet/
Reference:   http://russellgoyder.ca/geometric-algebra-cheat-sheet/

Covers sections 1-7 of the cheat sheet:
  1. Geometric product decomposition
  2. Vector-multivector products
  3. Vectors and bivectors
  4. Commutator product
  5. Pseudoscalar
  6. Outermorphism (rotation, determinant)
  7. Adjoint
"""
from sympy import symbols, cos, sin, S, simplify, trigsimp
from galgebra.ga import Ga

# 3D Euclidean GA: e1, e2, e3 with g = diag(1,1,1)
ga, e1, e2, e3 = Ga.build('e*1|2|3', g=[1, 1, 1])
I = ga.I()

a1, a2, a3 = symbols('a1 a2 a3', real=True)
b1, b2, b3 = symbols('b1 b2 b3', real=True)
c1, c2, c3 = symbols('c1 c2 c3', real=True)
d1, d2, d3 = symbols('d1 d2 d3', real=True)

a = a1*e1 + a2*e2 + a3*e3
b = b1*e1 + b2*e2 + b3*e3
c = c1*e1 + c2*e2 + c3*e3
d = d1*e1 + d2*e2 + d3*e3

print("=== Section 1: Geometric product decomposition ===")
# ab = a·b + a∧b
lhs = a * b
rhs = (a | b) + (a ^ b)
assert simplify((lhs - rhs).obj) == 0
print("ab = a·b + a∧b: OK")

# a·b = (ab + ba)/2  (symmetric part)
assert simplify(((a | b) - ga.mv(S.Half) * (a*b + b*a)).obj) == 0
print("a·b = (ab + ba)/2: OK")

# a∧b = (ab - ba)/2  (antisymmetric part)
assert simplify(((a ^ b) - ga.mv(S.Half) * (a*b - b*a)).obj) == 0
print("a∧b = (ab - ba)/2: OK")

print()
print("=== Section 2: Vector-multivector products ===")
# a*B_r = a<B_r + a∧B_r  (left contraction + outer product)
# For B_r = bivector (b^c):
B = b ^ c
lhs = a * B
rhs = (a < B) + (a ^ B)
assert simplify((lhs - rhs).obj) == 0
print("a*(b∧c) = a<(b∧c) + a∧b∧c: OK")

# For B_r = trivector (b^c^d):
T = b ^ c ^ d
lhs = a * T
rhs = (a < T) + (a ^ T)
assert simplify((lhs - rhs).obj) == 0
print("a*(b∧c∧d) = a<(b∧c∧d) + a∧b∧c∧d: OK")

print()
print("=== Section 3: Vectors and bivectors ===")
# a·(b∧c) = (a·b)c - (a·c)b
lhs = a < (b ^ c)
rhs = (a | b)*c - (a | c)*b
assert simplify((lhs - rhs).obj) == 0
print("a·(b∧c) = (a·b)c - (a·c)b: OK")

# (a∧b)·B = a·(b·B)  — left contraction associates
lhs = (a ^ b) < (c ^ d)
rhs = a < (b < (c ^ d))
assert simplify((lhs - rhs).obj) == 0
print("(a∧b)·B = a·(b·B): OK")

# (a∧b)·(c∧d) = (b·c)(a·d) - (a·c)(b·d)
lhs = (a ^ b) | (c ^ d)
rhs = (b | c)*(a | d) - (a | c)*(b | d)
assert simplify((lhs - rhs).obj) == 0
print("(a∧b)·(c∧d) = (b·c)(a·d) - (a·c)(b·d): OK")

print()
print("=== Section 4: Commutator product ===")
# M×N = (MN - NM)/2
# For bivector B and vector v: B×v = B·v
B2 = e1 ^ e2
for v in [e1, e2, e3, a]:
    comm = ga.mv(S.Half) * (B2*v - v*B2)
    inner = B2 | v
    assert simplify((comm - inner).obj) == 0
print("B×v = B·v (for all test vectors): OK")

# Commutator is grade-preserving for bivectors with bivectors
B3 = e2 ^ e3
comm_BB = ga.mv(S.Half) * (B2*B3 - B3*B2)
# result should be a bivector (grade 2)
assert comm_BB.pure_grade() == 2
print("B1×B2 is a bivector (grade-preserving): OK")

# Jacobi identity: A×(B×C) + B×(C×A) + C×(A×C) = 0  (using basis bivectors)
def comm(X, Y):
    return ga.mv(S.Half) * (X*Y - Y*X)

B4 = e1 ^ e3
jacobi = comm(B2, comm(B3, B4)) + comm(B3, comm(B4, B2)) + comm(B4, comm(B2, B3))
assert simplify(jacobi.obj) == 0
print("Jacobi identity for bivectors: OK")

print()
print("=== Section 5: Pseudoscalar ===")
# I^2 = -1 in 3D Euclidean
assert (I*I).scalar() == S.NegativeOne
print("I^2 = -1: OK")

# IA_r = (-1)^(r(n-1)) * A_r*I  (n=3, so (-1)^(2r) = 1 for all r)
# I commutes with all grades in 3D Euclidean
for mv, name in [(e1, "e1"), (e1^e2, "e1^e2"), (I, "I")]:
    assert simplify((I*mv - mv*I).obj) == 0
print("I commutes with all blades in 3D ((-1)^(2r)=1): OK")

# dual via I: dual(A) = A * I^{-1} = -A*I (since I^{-1} = -I when I^2=-1)
assert simplify((I.inv() + I).obj) == 0, "I^{-1} = -I when I^2=-1"
e1_dual = e1 * I.inv()
assert e1_dual == -(e2 ^ e3)
print("dual(e1) = e1*I^{-1} = -e2^e3: OK")

print()
print("=== Section 6: Outermorphism / Rotation ===")
theta = symbols('theta', real=True)
B_hat = e1 ^ e2  # unit bivector (already unit in Euclidean)

# Rotor: R = exp(-B*theta/2) = cos(theta/2) - sin(theta/2)*B
R = ga.mv(cos(theta/2)) - ga.mv(sin(theta/2)) * B_hat
R_rev = ~R

# R*~R = 1
assert trigsimp((R * R_rev).scalar()) == 1
print("R*~R = 1: OK")

# Rotation of e1 by theta in e1^e2 plane
rotated = R * e1 * R_rev
rotated = ga.mv(trigsimp(rotated.obj))
expected = ga.mv(cos(theta))*e1 + ga.mv(sin(theta))*e2
assert trigsimp((rotated - expected).obj) == 0
print("R*e1*~R = cos(theta)*e1 + sin(theta)*e2: OK")

# Determinant: L(I) = det(L)*I  i.e. det = (L(I)*I^{-1}).scalar()
L = ga.lt(lambda x: R*x*R_rev)
det_val = trigsimp(L.det())
assert det_val == 1  # rotations have det=1
print("det(rotation) = 1: OK")

print()
print("=== Section 7: Adjoint ===")
# Adjoint defined by: a·F̄(b) = F(a)·b  i.e.  a | L(b) = L.adj()(a) | b
# Use a shear transformation for a non-trivial test
L_shear = ga.lt([e1, e1 + e2, e3])
L_adj = L_shear.adj()
lhs = a | L_shear(b)
rhs = L_adj(a) | b
assert simplify((lhs - rhs).obj) == 0
print("a·F̄(b) = F(a)·b: OK")

# det(FG) = det(F)*det(G)
phi = symbols('phi', real=True)
R2 = ga.mv(cos(phi/2)) - ga.mv(sin(phi/2)) * (e2 ^ e3)
L2 = ga.lt(lambda x: R2*x*~R2)

def lt_compose(f, g):
    return ga.lt(lambda x: f(g(x)))

L_composed = lt_compose(L, L2)
det_product = trigsimp(L.det() * L2.det())
det_composed = trigsimp(L_composed.det())
assert trigsimp(det_product - det_composed) == 0
print("det(FG) = det(F)*det(G): OK")

print()
print("All cheat sheet identities validated successfully!")
