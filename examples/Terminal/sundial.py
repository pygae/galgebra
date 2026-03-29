"""
Sundial problem: core geometric algebra ideas.

Inspired by Russell Goyder's sundial analysis using GA:
  - Paper: https://russellgoyder.github.io/sundial-latex/
  - Notebook: https://russellgoyder.github.io/sundial/
  - Cheat sheet: https://russellgoyder.github.io/geometric-algebra-cheat-sheet/

This script validates galgebra against the core GA operations used in the
sundial problem: rotations via rotors, reflections, and projections in
3D Euclidean geometric algebra.

Mirrors issue 506 from upstream.
"""
from sympy import symbols, cos, sin, pi, simplify, trigsimp, S
from galgebra.ga import Ga

# Set up 3D Euclidean GA
ga, e1, e2, e3 = Ga.build('e*1|2|3', g=[1, 1, 1])
I = ga.I()

# === Rotor: rotation in the e1^e2 plane by angle theta ===
theta = symbols('theta', real=True)
# Rotor R = exp(-theta/2 * e12) = cos(theta/2) - sin(theta/2) * e12
e12 = e1 ^ e2
R = ga.mv(cos(theta / 2)) - ga.mv(sin(theta / 2)) * e12
R_rev = ~R  # reverse

# Verify R * ~R = 1 (rotor normalization)
product = R * R_rev
assert trigsimp(product.scalar()) == 1, f"R*~R should be 1, got {product}"

# Apply rotation to e1: should give cos(theta)*e1 + sin(theta)*e2
rotated = R * e1 * R_rev
rotated_simplified = ga.mv(trigsimp(rotated.obj))
expected = ga.mv(cos(theta)) * e1 + ga.mv(sin(theta)) * e2
assert trigsimp((rotated_simplified - expected).obj) == 0, \
    f"Rotation failed: got {rotated_simplified}, expected {expected}"
print("Rotation of e1 by theta in e12 plane:", rotated_simplified)

# === Reflection: reflect vector through a plane ===
# Reflection of v through plane perpendicular to n: -n*v*n^{-1}
# For unit n: -n*v*n
n = e3  # reflect through xy-plane (perpendicular to e3)
a, b, c = symbols('a b c', real=True)
v = a * e1 + b * e2 + c * e3
reflected = -n * v * n
print("\nReflect (a*e1 + b*e2 + c*e3) through xy-plane:", reflected)
# Should negate the e3 component
expected_refl = a * e1 + b * e2 - c * e3
assert simplify((reflected - expected_refl).obj) == 0

# === Projection of vector onto a blade ===
# proj_B(v) = (v . B) * B^{-1}
# Project v onto the e1^e2 plane
B = e1 ^ e2
v_full = e1 + 2 * e2 + 3 * e3
# Using the contraction: (v < B) * B^{-1}
proj = (v_full < B) * B.inv()
print("\nProject (e1 + 2*e2 + 3*e3) onto e12 plane:", proj)
assert proj == e1 + 2 * e2

# === Dual: convert between blades and their complements ===
# In 3D, dual of e1 is e2^e3 (up to sign depending on convention)
# dual(e1) = e1 * I^{-1} = e1 * (-I) = -e1*I
# In 3D: I = e123, I^{-1} = -e123
# e1 * (-e123) = -e1*e1*e2*e3 = -e2*e3  (since e1*e1=1)
e1_dual = e1 * I.inv()
print("\nDual of e1:", e1_dual)
# dual(e1) = -e23 in the convention dual(A) = A * I^{-1}
assert e1_dual == -(e2 ^ e3)

# === Composition of rotations ===
phi = symbols('phi', real=True)
R1 = ga.mv(cos(theta / 2)) - ga.mv(sin(theta / 2)) * (e1 ^ e2)
R2 = ga.mv(cos(phi / 2)) - ga.mv(sin(phi / 2)) * (e2 ^ e3)

# The composed rotation R2*R1 is also a rotor
R_composed = R2 * R1
R_composed_rev = ~R_composed
norm_sq = trigsimp((R_composed * R_composed_rev).scalar())
assert norm_sq == 1, f"Composed rotor not normalized: {norm_sq}"
print("\nComposed rotor R2*R1 is normalized: True")

# === Sundial core idea: hour angle to shadow direction ===
# The gnomon direction is along e3 (pointing to celestial pole)
# The sun's hour angle h rotates the sun direction in the e1^e2 plane
# The shadow is the projection of the gnomon onto the dial plane
h = symbols('h', real=True)
R_sun = ga.mv(cos(h / 2)) - ga.mv(sin(h / 2)) * (e1 ^ e2)
sun_dir = R_sun * e1 * ~R_sun
sun_dir = ga.mv(trigsimp(sun_dir.obj))
print("\nSun direction at hour angle h:", sun_dir)

print("\nAll sundial core GA operations validated successfully!")
