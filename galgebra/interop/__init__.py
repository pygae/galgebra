"""
Interoperability interfaces for creating geometric algebras using
conventions from other libraries.

The ``Cl(p, q, r)`` signature interface was popularized by ganja.js and
adopted by kingdon.  Two flavours are provided:

- ``galgebra.interop.Cl`` uses galgebra defaults (no surprises for
  existing users).
- ``galgebra.interop.kingdon.Cl`` uses kingdon conventions
  (``dual_mode='Iinv+'``, 0-indexed basis names).

Known incompatibilities between galgebra and kingdon:

- **Basis naming**: galgebra numbers from 1 (``e_1, e_2, ...``);
  kingdon defaults to 0-indexed for PGA (``e0, e1, e2, e3``).
- **Dual convention**: galgebra's default is ``'I+'`` (``dual(A) = I * A``);
  kingdon uses ``'Iinv+'`` (``dual(A) = A * I^{-1}``) via its
  ``codegen_polarity``.
"""

from ._cl import Cl  # noqa: F401
