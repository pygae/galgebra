"""
Regression tests validating galgebra against Russell Goyder's geometric
algebra cheat sheet (https://russellgoyder.github.io/geometric-algebra-cheat-sheet/).

Covers all 7 sections in G(3,0): geometric product, vector-multivector
products, bivector identities, commutator product, pseudoscalar,
outermorphism/rotation, and adjoint.
"""
import pytest
from sympy import symbols, cos, sin, S, simplify, trigsimp, Symbol, Matrix

from galgebra.ga import Ga


@pytest.fixture(scope='module')
def g3():
    ga, e1, e2, e3 = Ga.build('e*1|2|3', g=[1, 1, 1])
    return ga, e1, e2, e3


@pytest.fixture(scope='module')
def vecs(g3):
    ga, e1, e2, e3 = g3
    a1, a2, a3 = symbols('a1 a2 a3', real=True)
    b1, b2, b3 = symbols('b1 b2 b3', real=True)
    c1, c2, c3 = symbols('c1 c2 c3', real=True)
    d1, d2, d3 = symbols('d1 d2 d3', real=True)
    a = a1*e1 + a2*e2 + a3*e3
    b = b1*e1 + b2*e2 + b3*e3
    c = c1*e1 + c2*e2 + c3*e3
    d = d1*e1 + d2*e2 + d3*e3
    return a, b, c, d


def test_geometric_product_decomposition(g3, vecs):
    """Section 1: ab = a·b + a∧b, inner/outer as symmetric/antisymmetric parts."""
    ga, e1, e2, e3 = g3
    a, b, c, d = vecs

    assert simplify((a*b - (a|b) - (a^b)).obj) == 0
    assert simplify(((a|b) - ga.mv(S.Half)*(a*b + b*a)).obj) == 0
    assert simplify(((a^b) - ga.mv(S.Half)*(a*b - b*a)).obj) == 0


def test_vector_multivector_products(g3, vecs):
    """Section 2: a*B = a<B + a∧B (grade raising/lowering)."""
    ga, e1, e2, e3 = g3
    a, b, c, d = vecs

    B = b ^ c
    assert simplify((a*B - (a < B) - (a ^ B)).obj) == 0

    T = b ^ c ^ d
    assert simplify((a*T - (a < T) - (a ^ T)).obj) == 0


def test_bivector_identities(g3, vecs):
    """Section 3: a·(b∧c), (a∧b)·B associativity, (a∧b)·(c∧d) scalar formula."""
    ga, e1, e2, e3 = g3
    a, b, c, d = vecs

    # a·(b∧c) = (a·b)c - (a·c)b
    lhs1 = a < (b^c)
    rhs1 = (a|b)*c - (a|c)*b
    assert simplify((lhs1 - rhs1).obj) == 0

    # (a∧b)·B = a·(b·B)  (left contraction associates)
    lhs2 = (a^b) < (c^d)
    rhs2 = a < (b < (c^d))
    assert simplify((lhs2 - rhs2).obj) == 0

    # (a∧b)·(c∧d) = (b·c)(a·d) - (a·c)(b·d)
    lhs3 = (a^b) | (c^d)
    rhs3 = (b|c)*(a|d) - (a|c)*(b|d)
    assert simplify((lhs3 - rhs3).obj) == 0


def test_commutator_product(g3, vecs):
    """Section 4: B×v = B·v, grade-preserving on bivectors, Jacobi identity."""
    ga, e1, e2, e3 = g3
    a, b, c, d = vecs

    def comm(X, Y):
        return ga.mv(S.Half) * (X*Y - Y*X)

    B2 = e1 ^ e2
    for v in [e1, e2, e3, a]:
        assert simplify((comm(B2, v) - (B2 | v)).obj) == 0

    B3 = e2 ^ e3
    assert (comm(B2, B3)).pure_grade() == 2

    B4 = e1 ^ e3
    jacobi = comm(B2, comm(B3, B4)) + comm(B3, comm(B4, B2)) + comm(B4, comm(B2, B3))
    assert simplify(jacobi.obj) == 0


def test_pseudoscalar(g3):
    """Section 5: I²=-1, I central in G(3,0), dual via I⁻¹."""
    ga, e1, e2, e3 = g3
    I = ga.I()

    assert (I*I).scalar() == S.NegativeOne

    # (-1)^(r*(n-1)) = (-1)^(2r) = 1 for all r in 3D: I is central
    for mv in [e1, e1^e2, I]:
        assert simplify((I*mv - mv*I).obj) == 0

    # dual(A) = A * I^{-1}  (Iinv+ convention, not galgebra's built-in dual()
    # which uses I+ i.e. dual(A) = I*A — a different sign convention)
    assert simplify((I.inv() + I).obj) == 0   # I^{-1} = -I when I^2=-1
    assert e1 * I.inv() == -(e2 ^ e3)


def test_outermorphism_rotation(g3):
    """Section 6: rotor normalization, rotation sandwich, det(rotation)=1."""
    ga, e1, e2, e3 = g3
    theta = symbols('theta', real=True)
    R = ga.mv(cos(theta/2)) - ga.mv(sin(theta/2)) * (e1^e2)

    assert trigsimp((R * ~R).scalar()) == 1

    rotated = ga.mv(trigsimp((R * e1 * ~R).obj))
    assert trigsimp((rotated - ga.mv(cos(theta))*e1 - ga.mv(sin(theta))*e2).obj) == 0

    L = ga.lt(lambda x: R*x*~R)
    assert trigsimp(L.det()) == 1


def test_adjoint(g3, vecs):
    """Section 7: a·F̄(b) = F(a)·b and det(FG) = det(F)·det(G)."""
    ga, e1, e2, e3 = g3
    a, b, c, d = vecs

    L = ga.lt([e1, e1 + e2, e3])
    assert simplify((a | L(b) - L.adj()(a) | b).obj) == 0

    theta, phi = symbols('theta phi', real=True)
    R1 = ga.mv(cos(theta/2)) - ga.mv(sin(theta/2)) * (e1^e2)
    R2 = ga.mv(cos(phi/2)) - ga.mv(sin(phi/2)) * (e2^e3)
    L1 = ga.lt(lambda x: R1*x*~R1)
    L2 = ga.lt(lambda x: R2*x*~R2)
    L12 = ga.lt(lambda x: L1(L2(x)))
    assert trigsimp(L12.det() - L1.det()*L2.det()) == 0
